/*
 * The MIT License (MIT)
 * Copyright (c) 2023 Omics Data Automation, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.genomicsdb.importer.extensions;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.Message;
import com.google.protobuf.util.JsonFormat;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.genomicsdb.GenomicsDBUtils;
import org.genomicsdb.exception.GenomicsDBException;
import org.genomicsdb.model.Coordinates;
import org.genomicsdb.model.GenomicsDBCallsetsMapProto;
import org.genomicsdb.model.GenomicsDBImportConfiguration;
import org.genomicsdb.model.GenomicsDBVidMapProto;

import java.io.*;
import java.util.Base64;
import java.util.List;
import java.util.Set;

import static java.util.stream.Collectors.toList;

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
      String loaderJSONString = JsonFormat.printer().preservingProtoFieldNames().print(importConfiguration);

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
            throws GenomicsDBException {
      String callsetMapJSONString = null;
      try {
        callsetMapJSONString = JsonFormat.printer().preservingProtoFieldNames().print(callsetMappingPB);
      } catch (InvalidProtocolBufferException e) {
        throw new GenomicsDBException(String.format("Exception while printing callsetMapping protobuf message: %s", e.getLocalizedMessage()));
      }
      if (GenomicsDBUtils.writeToFile(outputCallsetMapJSONFilePath, callsetMapJSONString) != 0) {
        throw new GenomicsDBException(String.format("Could not write callset map json file : %s", outputCallsetMapJSONFilePath));
      }
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
            throws GenomicsDBException {
      String vidMapJSONString;
      try {
        vidMapJSONString = JsonFormat.printer().preservingProtoFieldNames().print(vidMappingPB);
      } catch (InvalidProtocolBufferException e) {
        throw new GenomicsDBException(String.format("Could not get string form for vidMappingPB: %s", e.getLocalizedMessage()));
      }
      if (GenomicsDBUtils.writeToFile(outputVidMapJSONFilePath, vidMapJSONString) != 0) {
        throw new GenomicsDBException(String.format("Could not write vid map json file : %s", outputVidMapJSONFilePath));
      }
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
            throws GenomicsDBException {
      final OutputStream stream = new ByteArrayOutputStream();
      final VCFHeader vcfHeader = new VCFHeader(headerLines);
      VariantContextWriter vcfWriter = new VariantContextWriterBuilder()
              .clearOptions()
              .setOutputStream(stream)
              .build();
      vcfWriter.writeHeader(vcfHeader);
      vcfWriter.close();
      String buffer = stream.toString();
      if (GenomicsDBUtils.writeToFile(outputVcfHeaderFilePath, buffer) != 0) {
        throw new GenomicsDBException(String.format("Could not write output header file : %s", outputVcfHeaderFilePath));
      }
    }

    /**
     * Generate the ProtoBuf data structure for vid mapping
     * from the existing vid file
     *
     * @param vidFile file with existing vid info
     * @return a vid map containing all field names, lengths and types
     * from the merged GVCF header. for incremental import case
     */
    default GenomicsDBVidMapProto.VidMappingPB generateVidMapFromFile(final String vidFile)
           throws GenomicsDBException {
      String existingVidJson = GenomicsDBUtils.readEntireFile(vidFile);
      GenomicsDBVidMapProto.VidMappingPB.Builder vidMapBuilder =
              GenomicsDBVidMapProto.VidMappingPB.newBuilder();
      try {
        JsonFormat.parser().merge(existingVidJson, vidMapBuilder);
      } catch (IOException e) {
        throw new GenomicsDBException(String.format("Could not generate vidmap from file : %s", e.getMessage()));
      }
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
    default boolean checkVidIfContigStartsWithinColumnBounds(GenomicsDBVidMapProto.VidMappingPB vidmap,
                                                             long[] bounds, Coordinates.ContigInterval contigInterval)
            throws GenomicsDBException {
        List<GenomicsDBVidMapProto.Chromosome> contigList = vidmap.getContigsList().stream().filter(
            x -> x.getName().equals(contigInterval.getContig())).collect(toList());
        if (contigList.size() != 1) {
            throw new GenomicsDBException("Contig "+contigInterval.getContig()+" not found, or found multiple times in vid");
        }
        return contigList.get(0).getTiledbColumnOffset() < bounds[1] &&
            contigList.get(0).getTiledbColumnOffset() >= bounds[0];
    }

    /**
     * Serialize Protobuf json file into a base64 encoded string
     * 
     * @param builder protobuf message builder
     * @param file Protobuf jsonfile
     * @return String base64 encoded protobuf string
     */
    static String getProtobufAsBase64StringFromFile(Message.Builder builder,
            String file) throws GenomicsDBException{
        String jsonString = GenomicsDBUtils.readEntireFile(file);
        try {
          JsonFormat.parser().ignoringUnknownFields().merge(jsonString, builder);
        } catch (IOException e) {
          throw new GenomicsDBException(String.format("Serialize protobuf json file into base64 encoded encoded string: %s", e.getMessage()));
        }
        byte[] pb = builder.build().toByteArray();
        return Base64.getEncoder().encodeToString(pb);
    }

    /**
     * Create Protobuf message from base64 encoded string
     * @param builder protobuf message builder
     * @param pbString base64 encoded protobug string
     * @return protobuf message
     * @throws InvalidProtocolBufferException when protobuf cannot parse the string
     */
    static Message getProtobufFromBase64EncodedString(Message.Builder builder, String pbString) 
            throws InvalidProtocolBufferException {
        byte[] pbDecoded = Base64.getDecoder().decode(pbString);
        builder.mergeFrom(pbDecoded);
        return builder.build();
    }
}
