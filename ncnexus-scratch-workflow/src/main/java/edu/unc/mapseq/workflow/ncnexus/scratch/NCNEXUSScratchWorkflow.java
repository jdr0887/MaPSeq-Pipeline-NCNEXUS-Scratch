package edu.unc.mapseq.workflow.ncnexus.scratch;

import java.io.File;
import java.util.List;
import java.util.ResourceBundle;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobBuilder;
import org.renci.jlrm.condor.CondorJobEdge;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.model.Flowcell;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.WorkflowRunAttempt;
import edu.unc.mapseq.module.bwa.BWAAlignCLI;
import edu.unc.mapseq.module.bwa.BWASAMPairedEndCLI;
import edu.unc.mapseq.module.core.RemoveCLI;
import edu.unc.mapseq.module.core.WriteVCFHeaderCLI;
import edu.unc.mapseq.module.fastqc.FastQCCLI;
import edu.unc.mapseq.module.fastqc.IgnoreLevelType;
import edu.unc.mapseq.module.gatk.GATKDownsamplingType;
import edu.unc.mapseq.module.gatk.GATKPhoneHomeType;
import edu.unc.mapseq.module.gatk2.GATKDepthOfCoverageCLI;
import edu.unc.mapseq.module.gatk2.GATKUnifiedGenotyperCLI;
import edu.unc.mapseq.module.picard.PicardAddOrReplaceReadGroupsCLI;
import edu.unc.mapseq.module.picard.PicardMarkDuplicatesCLI;
import edu.unc.mapseq.module.picard.PicardSortOrderType;
import edu.unc.mapseq.module.samtools.SAMToolsFlagstatCLI;
import edu.unc.mapseq.module.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.impl.AbstractSampleWorkflow;
import edu.unc.mapseq.workflow.impl.WorkflowJobFactory;
import edu.unc.mapseq.workflow.impl.WorkflowUtil;

public class NCNEXUSScratchWorkflow extends AbstractSampleWorkflow {

    private final Logger logger = LoggerFactory.getLogger(NCNEXUSScratchWorkflow.class);

    public NCNEXUSScratchWorkflow() {
        super();
    }

    @Override
    public String getName() {
        return NCNEXUSScratchWorkflow.class.getSimpleName().replace("Workflow", "");
    }

    @Override
    public String getVersion() {
        ResourceBundle bundle = ResourceBundle.getBundle("edu/unc/mapseq/workflow/ncnexus/scratch/workflow");
        String version = bundle.getString("version");
        return StringUtils.isNotEmpty(version) ? version : "0.0.1-SNAPSHOT";
    }

    @Override
    public Graph<CondorJob, CondorJobEdge> createGraph() throws WorkflowException {
        logger.info("ENTERING createGraph()");

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

        Set<Sample> sampleSet = getAggregatedSamples();
        logger.info("sampleSet.size(): {}", sampleSet.size());

        String siteName = getWorkflowBeanService().getAttributes().get("siteName");
        String referenceSequence = getWorkflowBeanService().getAttributes().get("referenceSequence");
        String readGroupPlatform = getWorkflowBeanService().getAttributes().get("readGroupPlatform");
        String readGroupPlatformUnit = getWorkflowBeanService().getAttributes().get("readGroupPlatformUnit");
        String depthOfCoverageIntervalList = getWorkflowBeanService().getAttributes()
                .get("depthOfCoverageIntervalList");
        String unifiedGenotyperIntervalList = getWorkflowBeanService().getAttributes().get(
                "unifiedGenotyperIntervalList");
        String unifiedGenotyperDBSNP = getWorkflowBeanService().getAttributes().get("unifiedGenotyperDBSNP");
        String GATKKey = getWorkflowBeanService().getAttributes().get("GATKKey");

        WorkflowRunAttempt attempt = getWorkflowRunAttempt();

        for (Sample sample : sampleSet) {

            if ("Undetermined".equals(sample.getBarcode())) {
                continue;
            }

            logger.debug(sample.toString());

            Flowcell flowcell = sample.getFlowcell();
            File outputDirectory = new File(sample.getOutputDirectory(), getName());
            File tmpDirectory = new File(outputDirectory, "tmp");
            tmpDirectory.mkdirs();

            List<File> readPairList = WorkflowUtil.getReadPairList(sample.getFileDatas(), flowcell.getName(),
                    sample.getLaneIndex());
            logger.debug("readPairList.size(): {}", readPairList.size());

            // assumption: a dash is used as a delimiter between a participantId
            // and the external code
            int idx = sample.getName().lastIndexOf("-");
            String sampleName = idx != -1 ? sample.getName().substring(0, idx) : sample.getName();

            if (readPairList.size() != 2) {
                throw new WorkflowException("readPairList != 2");
            }

            File r1FastqFile = readPairList.get(0);
            String r1FastqRootName = WorkflowUtil.getRootFastqName(r1FastqFile.getName());

            File r2FastqFile = readPairList.get(1);
            String r2FastqRootName = WorkflowUtil.getRootFastqName(r2FastqFile.getName());

            String fastqLaneRootName = StringUtils.removeEnd(r2FastqRootName, "_R2");

            File writeVCFHeaderOut;
            File fastqcR1Output;
            File fastqcR2Output;
            try {
                // new job
                CondorJobBuilder builder = WorkflowJobFactory.createJob(++count, WriteVCFHeaderCLI.class,
                        attempt.getId(), sample.getId()).siteName(siteName);
                String flowcellProper = flowcell.getName().substring(flowcell.getName().length() - 9,
                        flowcell.getName().length());
                writeVCFHeaderOut = new File(outputDirectory, fastqLaneRootName + ".vcf.hdr");
                builder.addArgument(WriteVCFHeaderCLI.BARCODE, sample.getBarcode())
                        .addArgument(WriteVCFHeaderCLI.RUN, flowcell.getName())
                        .addArgument(WriteVCFHeaderCLI.SAMPLENAME, sampleName)
                        .addArgument(WriteVCFHeaderCLI.STUDYNAME, sample.getStudy().getName())
                        .addArgument(WriteVCFHeaderCLI.LANE, sample.getLaneIndex().toString())
                        .addArgument(WriteVCFHeaderCLI.LABNAME, "kwilhelmsen")
                        .addArgument(WriteVCFHeaderCLI.FLOWCELL, flowcellProper)
                        .addArgument(WriteVCFHeaderCLI.OUTPUT, writeVCFHeaderOut.getAbsolutePath());
                CondorJob writeVCFHeaderJob = builder.build();
                logger.info(writeVCFHeaderJob.toString());
                graph.addVertex(writeVCFHeaderJob);

                // new job
                builder = WorkflowJobFactory.createJob(++count, FastQCCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName);
                fastqcR1Output = new File(outputDirectory, r1FastqRootName + ".fastqc.zip");
                builder.addArgument(FastQCCLI.INPUT, r1FastqFile.getAbsolutePath())
                        .addArgument(FastQCCLI.OUTPUT, fastqcR1Output.getAbsolutePath())
                        .addArgument(FastQCCLI.IGNORE, IgnoreLevelType.ERROR.toString());
                CondorJob fastQCR1Job = builder.build();
                logger.info(fastQCR1Job.toString());
                graph.addVertex(fastQCR1Job);

                // new job
                builder = WorkflowJobFactory
                        .createJob(++count, BWAAlignCLI.class, attempt.getId(), sample.getId(), false)
                        .siteName(siteName).numberOfProcessors(4);
                File saiR1OutFile = new File(outputDirectory, r1FastqRootName + ".sai");
                builder.addArgument(BWAAlignCLI.THREADS, "4")
                        .addArgument(BWAAlignCLI.FASTQ, r1FastqFile.getAbsolutePath())
                        .addArgument(BWAAlignCLI.FASTADB, referenceSequence)
                        .addArgument(BWAAlignCLI.OUTFILE, saiR1OutFile.getAbsolutePath());
                CondorJob bwaAlignR1Job = builder.build();
                logger.info(bwaAlignR1Job.toString());
                graph.addVertex(bwaAlignR1Job);
                graph.addEdge(fastQCR1Job, bwaAlignR1Job);

                // new job
                builder = WorkflowJobFactory.createJob(++count, FastQCCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName);
                fastqcR2Output = new File(outputDirectory, r2FastqRootName + ".fastqc.zip");
                builder.addArgument(FastQCCLI.INPUT, r2FastqFile.getAbsolutePath())
                        .addArgument(FastQCCLI.OUTPUT, fastqcR2Output.getAbsolutePath())
                        .addArgument(FastQCCLI.IGNORE, IgnoreLevelType.ERROR.toString());
                CondorJob fastQCR2Job = builder.build();
                logger.info(fastQCR2Job.toString());
                graph.addVertex(fastQCR2Job);

                // new job
                builder = WorkflowJobFactory
                        .createJob(++count, BWAAlignCLI.class, attempt.getId(), sample.getId(), false)
                        .siteName(siteName).numberOfProcessors(4);
                File saiR2OutFile = new File(outputDirectory, r2FastqRootName + ".sai");
                builder.addArgument(BWAAlignCLI.THREADS, "4")
                        .addArgument(BWAAlignCLI.FASTQ, r2FastqFile.getAbsolutePath())
                        .addArgument(BWAAlignCLI.FASTADB, referenceSequence)
                        .addArgument(BWAAlignCLI.OUTFILE, saiR2OutFile.getAbsolutePath());
                CondorJob bwaAlignR2Job = builder.build();
                logger.info(bwaAlignR2Job.toString());
                graph.addVertex(bwaAlignR2Job);
                graph.addEdge(fastQCR2Job, bwaAlignR2Job);

                // new job
                builder = WorkflowJobFactory.createJob(++count, BWASAMPairedEndCLI.class, attempt.getId(),
                        sample.getId(), false).siteName(siteName);
                File bwaSAMPairedEndOutFile = new File(outputDirectory, fastqLaneRootName + ".sam");
                builder.addArgument(BWASAMPairedEndCLI.FASTADB, referenceSequence)
                        .addArgument(BWASAMPairedEndCLI.FASTQ1, r1FastqFile.getAbsolutePath())
                        .addArgument(BWASAMPairedEndCLI.FASTQ2, r2FastqFile.getAbsolutePath())
                        .addArgument(BWASAMPairedEndCLI.SAI1, saiR1OutFile.getAbsolutePath())
                        .addArgument(BWASAMPairedEndCLI.SAI2, saiR2OutFile.getAbsolutePath())
                        .addArgument(BWASAMPairedEndCLI.OUTFILE, bwaSAMPairedEndOutFile.getAbsolutePath());
                CondorJob bwaSAMPairedEndJob = builder.build();
                logger.info(bwaSAMPairedEndJob.toString());
                graph.addVertex(bwaSAMPairedEndJob);
                graph.addEdge(bwaAlignR1Job, bwaSAMPairedEndJob);
                graph.addEdge(bwaAlignR2Job, bwaSAMPairedEndJob);
                graph.addEdge(writeVCFHeaderJob, bwaSAMPairedEndJob);

                // new job
                builder = WorkflowJobFactory
                        .createJob(++count, RemoveCLI.class, attempt.getId(), sample.getId(), false).siteName(siteName);
                builder.addArgument(RemoveCLI.FILE, saiR1OutFile.getAbsolutePath()).addArgument(RemoveCLI.FILE,
                        saiR2OutFile.getAbsolutePath());
                CondorJob removeSAIJob = builder.build();
                logger.info(removeSAIJob.toString());
                graph.addVertex(removeSAIJob);
                graph.addEdge(bwaSAMPairedEndJob, removeSAIJob);

                // new job
                builder = WorkflowJobFactory.createJob(++count, PicardAddOrReplaceReadGroupsCLI.class, attempt.getId(),
                        sample.getId()).siteName(siteName);
                File picardAddOrReplaceReadGroupsOuput = new File(outputDirectory, bwaSAMPairedEndOutFile.getName()
                        .replace(".sam", ".fixed-rg.bam"));
                builder.addArgument(PicardAddOrReplaceReadGroupsCLI.INPUT, bwaSAMPairedEndOutFile.getAbsolutePath())
                        .addArgument(PicardAddOrReplaceReadGroupsCLI.OUTPUT,
                                picardAddOrReplaceReadGroupsOuput.getAbsolutePath())
                        .addArgument(PicardAddOrReplaceReadGroupsCLI.SORTORDER,
                                PicardSortOrderType.COORDINATE.toString().toLowerCase())
                        .addArgument(
                                PicardAddOrReplaceReadGroupsCLI.READGROUPID,
                                String.format("%s-%s_L%03d", flowcell.getName(), sample.getBarcode(),
                                        sample.getLaneIndex()))
                        .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPLIBRARY, sampleName)
                        .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORM, readGroupPlatform)
                        .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORMUNIT, sample.getBarcode())
                        .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPSAMPLENAME, sampleName)
                        .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPCENTERNAME, "UNC");
                CondorJob picardAddOrReplaceReadGroupsJob = builder.build();
                logger.info(picardAddOrReplaceReadGroupsJob.toString());
                graph.addVertex(picardAddOrReplaceReadGroupsJob);
                graph.addEdge(bwaSAMPairedEndJob, picardAddOrReplaceReadGroupsJob);

                // new job
                builder = WorkflowJobFactory
                        .createJob(++count, RemoveCLI.class, attempt.getId(), sample.getId(), false).siteName(siteName);
                builder.addArgument(RemoveCLI.FILE, bwaSAMPairedEndOutFile.getAbsolutePath());
                CondorJob removeBWASAMPairedEndOutFileJob = builder.build();
                logger.info(removeBWASAMPairedEndOutFileJob.toString());
                graph.addVertex(removeBWASAMPairedEndOutFileJob);
                graph.addEdge(picardAddOrReplaceReadGroupsJob, removeBWASAMPairedEndOutFileJob);

                // new job
                builder = WorkflowJobFactory
                        .createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
                File samtoolsIndexOutput = new File(outputDirectory, picardAddOrReplaceReadGroupsOuput.getName()
                        .replace(".bam", ".bai"));
                builder.addArgument(SAMToolsIndexCLI.INPUT, picardAddOrReplaceReadGroupsOuput.getAbsolutePath())
                        .addArgument(SAMToolsIndexCLI.OUTPUT, samtoolsIndexOutput.getAbsolutePath());
                CondorJob samtoolsIndexJob = builder.build();
                logger.info(samtoolsIndexJob.toString());
                graph.addVertex(samtoolsIndexJob);
                graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsIndexJob);

                // deduped job
                File dedupedBamFile = new File(outputDirectory, picardAddOrReplaceReadGroupsOuput.getName().replace(
                        ".bam", ".deduped.bam"));
                builder = WorkflowJobFactory.createJob(++count, PicardMarkDuplicatesCLI.class, attempt.getId(),
                        sample.getId()).siteName(siteName);
                File picardMarkDuplicatesMetricsFile = new File(outputDirectory, dedupedBamFile.getName().replace(
                        ".bam", ".metrics"));
                builder.addArgument(PicardMarkDuplicatesCLI.INPUT, picardAddOrReplaceReadGroupsOuput.getAbsolutePath())
                        .addArgument(PicardMarkDuplicatesCLI.OUTPUT, dedupedBamFile.getAbsolutePath())
                        .addArgument(PicardMarkDuplicatesCLI.METRICSFILE,
                                picardMarkDuplicatesMetricsFile.getAbsolutePath());
                CondorJob dedupedBamJob = builder.build();
                logger.info(dedupedBamJob.toString());
                graph.addVertex(dedupedBamJob);
                graph.addEdge(samtoolsIndexJob, dedupedBamJob);

                // index job
                builder = WorkflowJobFactory
                        .createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
                File dedupedBaiFile = new File(outputDirectory, dedupedBamFile.getName().replace(".bam", ".bai"));
                builder.addArgument(SAMToolsIndexCLI.INPUT, dedupedBamFile.getAbsolutePath()).addArgument(
                        SAMToolsIndexCLI.OUTPUT, dedupedBaiFile.getAbsolutePath());
                CondorJob dedupedBaiJob = builder.build();
                logger.info(dedupedBaiJob.toString());
                graph.addVertex(dedupedBaiJob);
                graph.addEdge(dedupedBamJob, dedupedBaiJob);

                // flagstat job
                builder = WorkflowJobFactory.createJob(++count, SAMToolsFlagstatCLI.class, attempt.getId(),
                        sample.getId()).siteName(siteName);
                File dedupedRealignFixPrintReadsFlagstatFile = new File(outputDirectory, dedupedBamFile.getName()
                        .replace(".bam", ".realign.fix.pr.flagstat"));
                builder.addArgument(SAMToolsFlagstatCLI.INPUT, dedupedBamFile.getAbsolutePath()).addArgument(
                        SAMToolsFlagstatCLI.OUTPUT, dedupedRealignFixPrintReadsFlagstatFile.getAbsolutePath());
                CondorJob dedupedRealignFixPrintReadsFlagstatJob = builder.build();
                logger.info(dedupedRealignFixPrintReadsFlagstatJob.toString());
                graph.addVertex(dedupedRealignFixPrintReadsFlagstatJob);
                graph.addEdge(dedupedBaiJob, dedupedRealignFixPrintReadsFlagstatJob);

                // depth of coverage job
                builder = WorkflowJobFactory
                        .createJob(++count, GATKDepthOfCoverageCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName).initialDirectory(outputDirectory.getAbsolutePath());
                builder.addArgument(GATKDepthOfCoverageCLI.INPUTFILE, dedupedBamFile.getAbsolutePath())
                        .addArgument(GATKDepthOfCoverageCLI.OUTPUTPREFIX,
                                dedupedBamFile.getName().replace(".bam", ".realign.fix.pr.coverage"))
                        .addArgument(GATKDepthOfCoverageCLI.KEY, GATKKey)
                        .addArgument(GATKDepthOfCoverageCLI.REFERENCESEQUENCE, referenceSequence)
                        .addArgument(GATKDepthOfCoverageCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                        .addArgument(GATKDepthOfCoverageCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                        .addArgument(GATKDepthOfCoverageCLI.VALIDATIONSTRICTNESS, "LENIENT")
                        .addArgument(GATKDepthOfCoverageCLI.OMITDEPTHOUTPUTATEACHBASE)
                        .addArgument(GATKDepthOfCoverageCLI.INTERVALS, depthOfCoverageIntervalList);
                CondorJob dedupedRealignFixPrintReadsCoverageJob = builder.build();
                graph.addVertex(dedupedRealignFixPrintReadsCoverageJob);
                graph.addEdge(dedupedBaiJob, dedupedRealignFixPrintReadsCoverageJob);

                // unified genotyper job

                builder = WorkflowJobFactory
                        .createJob(++count, GATKUnifiedGenotyperCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName).numberOfProcessors(4);
                File dedupedRealignFixPrintReadsVcfFile = new File(outputDirectory, dedupedBamFile.getName().replace(
                        ".bam", ".realign.fix.pr.vcf"));
                File gatkUnifiedGenotyperMetrics = new File(outputDirectory, dedupedBamFile.getName().replace(".bam",
                        ".metrics"));
                builder.addArgument(GATKUnifiedGenotyperCLI.INPUTFILE, dedupedBamFile.getAbsolutePath())
                        .addArgument(GATKUnifiedGenotyperCLI.OUT, dedupedRealignFixPrintReadsVcfFile.getAbsolutePath())
                        .addArgument(GATKUnifiedGenotyperCLI.KEY, GATKKey)
                        .addArgument(GATKUnifiedGenotyperCLI.INTERVALS, unifiedGenotyperIntervalList)
                        .addArgument(GATKUnifiedGenotyperCLI.REFERENCESEQUENCE, referenceSequence)
                        .addArgument(GATKUnifiedGenotyperCLI.DBSNP, unifiedGenotyperDBSNP)
                        .addArgument(GATKUnifiedGenotyperCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                        .addArgument(GATKUnifiedGenotyperCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                        .addArgument(GATKUnifiedGenotyperCLI.GENOTYPELIKELIHOODSMODEL, "BOTH")
                        .addArgument(GATKUnifiedGenotyperCLI.OUTPUTMODE, "EMIT_ALL_SITES")
                        .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "AlleleBalance")
                        .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "DepthOfCoverage")
                        .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "HomopolymerRun")
                        .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "MappingQualityZero")
                        .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "QualByDepth")
                        .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "RMSMappingQuality")
                        .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "HaplotypeScore")
                        .addArgument(GATKUnifiedGenotyperCLI.DOWNSAMPLETOCOVERAGE, "250")
                        .addArgument(GATKUnifiedGenotyperCLI.STANDCALLCONF, "4")
                        .addArgument(GATKUnifiedGenotyperCLI.STANDEMITCONF, "0")
                        .addArgument(GATKUnifiedGenotyperCLI.NUMTHREADS, "4")
                        .addArgument(GATKUnifiedGenotyperCLI.METRICS, gatkUnifiedGenotyperMetrics.getAbsolutePath());
                CondorJob dedupedRealignFixPrintReadsVcfJob = builder.build();
                graph.addVertex(dedupedRealignFixPrintReadsVcfJob);
                graph.addEdge(dedupedRealignFixPrintReadsCoverageJob, dedupedRealignFixPrintReadsVcfJob);

            } catch (Exception e) {
                throw new WorkflowException(e);
            }
        }

        return graph;
    }

}
