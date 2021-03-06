package edu.unc.mapseq.workflow;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import org.junit.Test;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobBuilder;
import org.renci.jlrm.condor.CondorJobEdge;
import org.renci.jlrm.condor.ext.CondorDOTExporter;
import org.renci.jlrm.condor.ext.CondorJobVertexNameProvider;

import edu.unc.mapseq.module.bwa.BWAAlignCLI;
import edu.unc.mapseq.module.bwa.BWASAMPairedEndCLI;
import edu.unc.mapseq.module.core.WriteVCFHeaderCLI;
import edu.unc.mapseq.module.fastqc.FastQCCLI;
import edu.unc.mapseq.module.picard.PicardAddOrReplaceReadGroupsCLI;
import edu.unc.mapseq.module.samtools.SAMToolsIndexCLI;

public class WorkflowTest {

    @Test
    public void createDot() {

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

        // new job
        CondorJob fastQCR1Job = new CondorJobBuilder().name(
                String.format("%s_%d", FastQCCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(fastQCR1Job);

        // new job
        CondorJob bwaAlignR1Job = new CondorJobBuilder().name(
                String.format("%s_%d", BWAAlignCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(bwaAlignR1Job);
        graph.addEdge(fastQCR1Job, bwaAlignR1Job);

        // new job
        CondorJob fastQCR2Job = new CondorJobBuilder().name(
                String.format("%s_%d", FastQCCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(fastQCR2Job);

        // new job
        CondorJob bwaAlignR2Job = new CondorJobBuilder().name(
                String.format("%s_%d", BWAAlignCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(bwaAlignR2Job);
        graph.addEdge(fastQCR2Job, bwaAlignR2Job);

        CondorJob writeVCFHeaderJob = new CondorJobBuilder().name(
                String.format("%s_%d", WriteVCFHeaderCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(writeVCFHeaderJob);

        // new job
        CondorJob bwaSAMPairedEndJob = new CondorJobBuilder().name(
                String.format("%s_%d", BWASAMPairedEndCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(bwaSAMPairedEndJob);
        graph.addEdge(bwaAlignR1Job, bwaSAMPairedEndJob);
        graph.addEdge(bwaAlignR2Job, bwaSAMPairedEndJob);
        graph.addEdge(writeVCFHeaderJob, bwaSAMPairedEndJob);

        // new job
        CondorJob picardAddOrReplaceReadGroupsJob = new CondorJobBuilder().name(
                String.format("%s_%d", PicardAddOrReplaceReadGroupsCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(picardAddOrReplaceReadGroupsJob);
        graph.addEdge(bwaSAMPairedEndJob, picardAddOrReplaceReadGroupsJob);

        // new job
        CondorJob samtoolsIndexJob = new CondorJobBuilder().name(
                String.format("%s_%d", SAMToolsIndexCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(samtoolsIndexJob);
        graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsIndexJob);

        CondorJobVertexNameProvider vnp = new CondorJobVertexNameProvider();
        CondorDOTExporter<CondorJob, CondorJobEdge> dotExporter = new CondorDOTExporter<CondorJob, CondorJobEdge>(vnp,
                vnp, null, null, null, null);
        File srcSiteResourcesImagesDir = new File("src/site/resources/images");
        if (!srcSiteResourcesImagesDir.exists()) {
            srcSiteResourcesImagesDir.mkdirs();
        }
        File dotFile = new File(srcSiteResourcesImagesDir, "workflow.dag.dot");
        try {
            FileWriter fw = new FileWriter(dotFile);
            dotExporter.export(fw, graph);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    @Test
    public void messingWithGraph() {

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

        // new job
        CondorJob fastQCR1Job = new CondorJobBuilder().name(
                String.format("%s_%d", FastQCCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(fastQCR1Job);

        // new job
        CondorJob bwaAlignR1Job = new CondorJobBuilder().name(
                String.format("%s_%d", BWAAlignCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(bwaAlignR1Job);
        graph.addEdge(fastQCR1Job, bwaAlignR1Job);

        // new job
        CondorJob fastQCR2Job = new CondorJobBuilder().name(
                String.format("%s_%d", FastQCCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(fastQCR2Job);

        // new job
        CondorJob bwaAlignR2Job = new CondorJobBuilder().name(
                String.format("%s_%d", BWAAlignCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(bwaAlignR2Job);
        graph.addEdge(fastQCR2Job, bwaAlignR2Job);

        CondorJob writeVCFHeaderJob = new CondorJobBuilder().name(
                String.format("%s_%d", WriteVCFHeaderCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(writeVCFHeaderJob);

        // new job
        CondorJob bwaSAMPairedEndJob = new CondorJobBuilder().name(
                String.format("%s_%d", BWASAMPairedEndCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(bwaSAMPairedEndJob);
        graph.addEdge(bwaAlignR1Job, bwaSAMPairedEndJob);
        graph.addEdge(bwaAlignR2Job, bwaSAMPairedEndJob);
        graph.addEdge(writeVCFHeaderJob, bwaSAMPairedEndJob);

        // new job
        CondorJob picardAddOrReplaceReadGroupsJob = new CondorJobBuilder().name(
                String.format("%s_%d", PicardAddOrReplaceReadGroupsCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(picardAddOrReplaceReadGroupsJob);
        graph.addEdge(bwaSAMPairedEndJob, picardAddOrReplaceReadGroupsJob);

        // new job
        CondorJob samtoolsIndexJob = new CondorJobBuilder().name(
                String.format("%s_%d", SAMToolsIndexCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(samtoolsIndexJob);
        graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsIndexJob);

        Iterator<CondorJob> iter = new DepthFirstIterator<CondorJob, CondorJobEdge>(graph);
        CondorJob vertex;
        while (iter.hasNext()) {
            vertex = iter.next();
            System.out.println("Vertex " + vertex.getName() + " is connected to: " + graph.edgesOf(vertex).toString());
        }

    }

}
