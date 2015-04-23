package edu.unc.mapseq.ws.ncnexus.scratch;

import java.util.List;

import javax.jws.WebMethod;
import javax.jws.WebParam;
import javax.jws.WebService;
import javax.jws.soap.SOAPBinding;
import javax.jws.soap.SOAPBinding.Style;
import javax.jws.soap.SOAPBinding.Use;
import javax.ws.rs.Consumes;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;
import javax.xml.ws.BindingType;

@BindingType(value = javax.xml.ws.soap.SOAPBinding.SOAP11HTTP_BINDING)
@WebService(targetNamespace = "http://scratch.ncnexus.ws.mapseq.unc.edu", serviceName = "NCNEXUSScratchService", portName = "NCNEXUSScratchPort")
@SOAPBinding(style = Style.DOCUMENT, use = Use.LITERAL, parameterStyle = SOAPBinding.ParameterStyle.WRAPPED)
@Path("/NCNEXUSScratchService/")
@Consumes(MediaType.APPLICATION_JSON)
@Produces(MediaType.APPLICATION_JSON)
public interface NCNEXUSScratchService {

    @GET
    @Path("/lookupQuantificationResults/{sampleId}")
    @WebMethod
    public QualityControlInfo lookupQuantificationResults(
            @PathParam("sampleId") @WebParam(name = "sampleId") Long sampleId);

    @GET
    @Path("/lookupQuantificationResultsByFlowcell/{flowcellId}")
    @WebMethod
    public List<QualityControlInfo> lookupQuantificationResultsByFlowcell(
            @PathParam("flowcellId") @WebParam(name = "flowcellId") Long flowcellId);

}
