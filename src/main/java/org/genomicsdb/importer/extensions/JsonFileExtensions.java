package org.genomicsdb.importer.extensions;

import com.googlecode.protobuf.format.JsonFormat;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import org.genomicsdb.GenomicsDBUtils;
import org.genomicsdb.model.Coordinates;
import org.genomicsdb.model.GenomicsDBImportConfiguration;
import org.genomicsdb.model.GenomicsDBCallsetsMapProto;
import org.genomicsdb.model.GenomicsDBVidMapProto;

import java.io.*;
import java.util.Set;
import java.util.List;

import static java.util.stream.Collectors.toList;
import static com.googlecode.protobuf.format.JsonFormat.*;
import static org.genomicsdb.GenomicsDBUtils.*;

public interface JsonFileExtensions {
    /**
     * Create a JSON file from the import configuration
     *
     * @param importConfiguration The configuration object
     * @param filename            File to dump the loader JSON to
     * @return New file (with the specified name) written to local storage
     * @throws FileNotFoundException PrintWriter throws this exception when file not found
     */
    default File dumpTemporaryLoaderJSONFile(final GenomicsDBImportConfiguration.ImportConfiguration importConfiguration,
                                             final String filename) throws IOException {
        String loaderJSONString = JsonFormat.printToString(importConfiguration);

        File tempLoaderJSONFile = (filename.isEmpty()) ?
                File.createTempFile("loader_", ".json") :
                new File(filename);

        if (filename.isEmpty()) tempLoaderJSONFile.deleteOnExit();

        PrintWriter out = new PrintWriter(tempLoaderJSONFile);
        out.println(loaderJSONString);
        out.close();
        return tempLoaderJSONFile;
    }

    /**
     * Writes a JSON file from a vidmap protobuf object. This
     * method is not called implicitly by GenomicsDB constructor.
     * It needs to be called explicitly if the user wants these objects
     * to be written. Called explicitly from GATK-4 GenomicsDBImport tool
     *
     * @param outputCallsetMapJSONFilePath Full path of file to be written
     * @param callsetMappingPB             Protobuf callset map object
     */
    default void writeCallsetMapJSONFile(final String outputCallsetMapJSONFilePath,
                                         final GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB)
    {
        String callsetMapJSONString = printToString(callsetMappingPB);
	GenomicsDBUtils.writeToFile(outputCallsetMapJSONFilePath, callsetMapJSONString);
    }

    /**
     * Writes a JSON file from a set of VCF header lines. This
     * method is not called implicitly by GenomicsDB constructor.
     * It needs to be called explicitly if the user wants these objects
     * to be written. Called explicitly from GATK-4 GenomicsDBImport tool
     *
     * @param outputVidMapJSONFilePath Full path of file to be written
     * @param vidMappingPB             ProtoBuf data structure for vid mapping
     */
    default void writeVidMapJSONFile(final String outputVidMapJSONFilePath,
                                     final GenomicsDBVidMapProto.VidMappingPB vidMappingPB)
    {
        String vidMapJSONString = printToString(vidMappingPB);
	GenomicsDBUtils.writeToFile(outputVidMapJSONFilePath, vidMapJSONString);
    }

    /**
     * Writes a VCF Header file from a set of VCF header lines. This
     * method is not called implicitly by GenomicsDB constructor.
     * It needs to be called explicitly if the user wants these objects
     * to be written. Called explicitly from GATK-4 GenomicsDBImport tool
     *
     * @param outputVcfHeaderFilePath Full path of file to be written
     * @param headerLines             Set of header lines
     */
    default void writeVcfHeaderFile(final String outputVcfHeaderFilePath, final Set<VCFHeaderLine> headerLines)
    {
        final OutputStream stream = new ByteArrayOutputStream();
	final VCFHeader vcfHeader = new VCFHeader(headerLines);
	VariantContextWriter vcfWriter = new VariantContextWriterBuilder()
	    .clearOptions()
	    .setOutputStream(stream)
	    .build();
	vcfWriter.writeHeader(vcfHeader);
	vcfWriter.close();
	String buffer = stream.toString();
        GenomicsDBUtils.writeToFile(outputVcfHeaderFilePath, buffer);
    }

    /**
     * Generate the ProtoBuf data structure for vid mapping
     * from the existing vid file
     *
     * @param vidFile file with existing vid info
     * @return a vid map containing all field names, lengths and types
     * from the merged GVCF header. for incremental import case
     * @throws ParseException when there is an error parsing existing vid json
     */
    default GenomicsDBVidMapProto.VidMappingPB generateVidMapFromFile(final String vidFile)
        throws com.googlecode.protobuf.format.JsonFormat.ParseException {
        String existingVidJson = GenomicsDBUtils.readEntireFile(vidFile);
        GenomicsDBVidMapProto.VidMappingPB.Builder vidMapBuilder = 
                GenomicsDBVidMapProto.VidMappingPB.newBuilder();
        merge(existingVidJson, vidMapBuilder);
        return vidMapBuilder.build();
    }

    /**
     * Check if contig starts within Tiledb column bounds
     *
     * @param vidmap Vid protobuf object with contig mapping information
     * @param bounds array with lower and upper column bounds
     * @param contigInterval contig interval
     * @return true if contig starts within column bounds
     */
    default boolean checkVidForContigColumnOffsetOverlap(GenomicsDBVidMapProto.VidMappingPB vidmap,
            long[] bounds, Coordinates.ContigInterval contigInterval) {
        List<GenomicsDBVidMapProto.Chromosome> contigList = vidmap.getContigsList().stream().filter(
            x -> x.getName().equals(contigInterval.getContig())).collect(toList());
        if (contigList.size() != 1) {
            throw new RuntimeException("Contig "+contigInterval.getContig()+" not found, or found multiple times in vid");
        }
        return contigList.get(0).getTiledbColumnOffset() < bounds[1] &&
            contigList.get(0).getTiledbColumnOffset() >= bounds[0];
    }
}
