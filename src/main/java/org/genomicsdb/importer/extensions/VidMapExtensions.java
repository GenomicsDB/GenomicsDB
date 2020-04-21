package org.genomicsdb.importer.extensions;

import htsjdk.variant.vcf.*;

import org.genomicsdb.importer.Constants;
import org.genomicsdb.model.GenomicsDBVidMapProto;
import org.genomicsdb.model.Coordinates;
import org.genomicsdb.GenomicsDBUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.Map;
import java.util.HashMap;

import static java.util.stream.Collectors.toList;
import static com.googlecode.protobuf.format.JsonFormat.*;
import static org.genomicsdb.GenomicsDBUtils.readEntireFile;

public interface VidMapExtensions {

    public static String[] vcfHeaderLineTypeString = new String[] {
      "FILTER", "INFO", "FORMAT" };
    public static int vcfHeaderLineFilterIdx = 0;
    public static int vcfHeaderLineInfoIdx = 1;
    public static int vcfHeaderLineFormatIdx = 2;

    public class DuplicateFieldInfoTrackClass {
      public String vcfType;
      public String vcfLengthDescriptor;
      public boolean isInitialized = false;
      public int idxInGenomicsDBFieldInfoList;
    }

    public class DuplicateVcfFieldNamesCheckResult {
	public boolean skipCreatingNewEntry = false;
	public String fieldName;

	public DuplicateVcfFieldNamesCheckResult(final boolean b, final String s) {
	    skipCreatingNewEntry = b;
	    fieldName = s;
	}
    }

    /**
     * Generate the ProtoBuf data structure for vid mapping
     * Also remember, contigs are 1-based which means
     * TileDB column offsets should start from 1
     *
     * @param mergedHeader Header from all input GVCFs are merged to
     *                     create one vid map for all
     * @return a vid map containing all field names, lengths and types
     * from the merged GVCF header
     */
    default GenomicsDBVidMapProto.VidMappingPB generateVidMapFromMergedHeader(final Set<VCFHeaderLine> mergedHeader) {
        List<GenomicsDBVidMapProto.GenomicsDBFieldInfo> infoFields = new ArrayList<>();
        List<GenomicsDBVidMapProto.Chromosome> contigs = new ArrayList<>();

        //ID field
        GenomicsDBVidMapProto.GenomicsDBFieldInfo.Builder IDFieldBuilder = GenomicsDBVidMapProto.GenomicsDBFieldInfo.newBuilder();
        GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.Builder lengthDescriptorComponentBuilder =
                GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.newBuilder();
        lengthDescriptorComponentBuilder.setVariableLengthDescriptor("var");
        IDFieldBuilder.setName("ID").addType("char").addLength(lengthDescriptorComponentBuilder.build());
        infoFields.add(IDFieldBuilder.build());

        long columnOffset = 0L;

	Map<String, List<DuplicateFieldInfoTrackClass>> fieldNameToDuplicateCheckInfo =
	    new HashMap<String,List<DuplicateFieldInfoTrackClass>>();

        for (VCFHeaderLine headerLine : mergedHeader) {
            GenomicsDBVidMapProto.GenomicsDBFieldInfo.Builder infoBuilder = GenomicsDBVidMapProto.GenomicsDBFieldInfo.newBuilder();
            if (headerLine instanceof VCFFormatHeaderLine) {
                VCFFormatHeaderLine formatHeaderLine = (VCFFormatHeaderLine) headerLine;

		final String formatFieldName = formatHeaderLine.getID();
                boolean isGT = formatHeaderLine.getID().equals(VCFConstants.GENOTYPE_KEY);
                String genomicsDBType = isGT ? "int" : formatHeaderLine.getType().toString();
                String genomicsDBLength = isGT ? "PP" : (formatHeaderLine.getType() == VCFHeaderLineType.String)
                        ? "var" : getVcfHeaderLength(formatHeaderLine);
                lengthDescriptorComponentBuilder.setVariableLengthDescriptor(genomicsDBLength);

		final DuplicateVcfFieldNamesCheckResult dupCheck = checkForDuplicateVcfFieldNames(
			formatFieldName,
			vcfHeaderLineFormatIdx,
			genomicsDBType, genomicsDBLength,
			infoFields,
			fieldNameToDuplicateCheckInfo
			);
		if(dupCheck.skipCreatingNewEntry)
		    continue;

                infoBuilder.setName(dupCheck.fieldName)
                        .addType(genomicsDBType)
                        .addLength(lengthDescriptorComponentBuilder.build());
		if(!formatFieldName.equals(dupCheck.fieldName)) //field renamed due to collisions
		    infoBuilder.setVcfName(formatFieldName);

		infoBuilder.addVcfFieldClass("FORMAT");

                GenomicsDBVidMapProto.GenomicsDBFieldInfo formatField = infoBuilder.build();
                infoFields.add(formatField);
            } else if (headerLine instanceof VCFInfoHeaderLine) {
                VCFInfoHeaderLine infoHeaderLine = (VCFInfoHeaderLine) headerLine;
                final String infoFieldName = infoHeaderLine.getID();
		final String genomicsDBLength = infoHeaderLine.getType() == VCFHeaderLineType.String ?
		    "var" : getVcfHeaderLength(infoHeaderLine);
		final DuplicateVcfFieldNamesCheckResult dupCheck = checkForDuplicateVcfFieldNames(
			infoFieldName,
			vcfHeaderLineInfoIdx,
			infoHeaderLine.getType().toString(),
			genomicsDBLength,
			infoFields,
			fieldNameToDuplicateCheckInfo
			);
		if(dupCheck.skipCreatingNewEntry)
		    continue;
                infoBuilder.setName(dupCheck.fieldName);
		if(!infoFieldName.equals(dupCheck.fieldName)) //field renamed due to collisions
		    infoBuilder.setVcfName(infoFieldName);
                //allele specific annotations
                if (Constants.R_LENGTH_HISTOGRAM_FIELDS_FLOAT_BINS.contains(dupCheck.fieldName)
                        || Constants.R_LENGTH_TWO_DIM_FLOAT_VECTOR_FIELDS.contains(dupCheck.fieldName)
                        || Constants.R_LENGTH_TWO_DIM_INT_VECTOR_FIELDS.contains(dupCheck.fieldName)
                        ) {
                    lengthDescriptorComponentBuilder.setVariableLengthDescriptor("R");
                    infoBuilder.addLength(lengthDescriptorComponentBuilder.build());
                    lengthDescriptorComponentBuilder.setVariableLengthDescriptor("var"); //ignored - can set anything here
                    infoBuilder.addLength(lengthDescriptorComponentBuilder.build());
                    infoBuilder.addVcfDelimiter("|");
                    infoBuilder.addVcfDelimiter(",");
                    if (Constants.R_LENGTH_HISTOGRAM_FIELDS_FLOAT_BINS.contains(dupCheck.fieldName)) {
                        //Each element of the vector is a tuple <float, int>
                        infoBuilder.addType("float");
                        infoBuilder.addType("int");
                        infoBuilder.setVCFFieldCombineOperation("histogram_sum");
                    } else {
                        infoBuilder.setVCFFieldCombineOperation("element_wise_sum");
                        if (Constants.R_LENGTH_TWO_DIM_FLOAT_VECTOR_FIELDS.contains(dupCheck.fieldName)) {
                            infoBuilder.addType("float");
                        } else if (Constants.R_LENGTH_TWO_DIM_INT_VECTOR_FIELDS.contains(dupCheck.fieldName)) {
                            infoBuilder.addType("int");
                        }
                    }
                } else {
                    lengthDescriptorComponentBuilder.setVariableLengthDescriptor(genomicsDBLength);
                    infoBuilder.addType(infoHeaderLine.getType().toString())
                            .addLength(lengthDescriptorComponentBuilder.build());
                }
		infoBuilder.addVcfFieldClass("INFO");
                GenomicsDBVidMapProto.GenomicsDBFieldInfo infoField = infoBuilder.build();
                infoFields.add(infoField);
            } else if (headerLine instanceof VCFFilterHeaderLine) {
                VCFFilterHeaderLine filterHeaderLine = (VCFFilterHeaderLine) headerLine;
		final String filterFieldName = filterHeaderLine.getID();
		final DuplicateVcfFieldNamesCheckResult dupCheck = checkForDuplicateVcfFieldNames(
			filterFieldName,
			vcfHeaderLineFilterIdx,
			"Integer",
			"1",
			infoFields,
			fieldNameToDuplicateCheckInfo
			);
		if(dupCheck.skipCreatingNewEntry)
		    continue;
                infoBuilder.setName(dupCheck.fieldName);
		if(!filterFieldName.equals(dupCheck.fieldName)) //field renamed due to collisions
		    infoBuilder.setVcfName(filterFieldName);
		lengthDescriptorComponentBuilder.setVariableLengthDescriptor("1");
		infoBuilder.addType("Integer")
                            .addLength(lengthDescriptorComponentBuilder.build());
                infoBuilder.addVcfFieldClass("FILTER");
                GenomicsDBVidMapProto.GenomicsDBFieldInfo filterField = infoBuilder.build();
                infoFields.add(filterField);
            } else if (headerLine instanceof VCFContigHeaderLine) {
                VCFContigHeaderLine contigHeaderLine = (VCFContigHeaderLine) headerLine;
                long length = contigHeaderLine.getSAMSequenceRecord().getSequenceLength();
                GenomicsDBVidMapProto.Chromosome.Builder contigBuilder = GenomicsDBVidMapProto.Chromosome.newBuilder();
                contigBuilder.setName(contigHeaderLine.getID()).setLength(length).setTiledbColumnOffset(columnOffset);
                columnOffset += length;
                GenomicsDBVidMapProto.Chromosome chromosome = contigBuilder.build();
                contigs.add(chromosome);
            }
        }

        GenomicsDBVidMapProto.VidMappingPB.Builder vidMapBuilder = GenomicsDBVidMapProto.VidMappingPB.newBuilder();

        return vidMapBuilder.addAllFields(infoFields).addAllContigs(contigs).build();
    }

    /**
     * Maps the "Number" from INFO or FORMAT fields in VCF header
     * to GenomicsDB lengths. If unbounded, the field becomes a
     * variable length attribute (discussed in more detail in TileDB
     * tutorial tiledb.org). If "A", "R" or "G", the field also
     * becomes a variable length attribute, but the order and length
     * depends on order/number of alleles or genotypes. Integers
     * remain integers
     *
     * @param headerLine Info or Format header line from VCF
     * @return VAR, A, R, G, or integer values from VCF
     */
    default String getVcfHeaderLength(final VCFHeaderLine headerLine) {
        VCFHeaderLineCount type;

        if (headerLine instanceof VCFFormatHeaderLine) {
            type = ((VCFFormatHeaderLine) headerLine).getCountType();
        } else {
            type = ((VCFInfoHeaderLine) headerLine).getCountType();
        }

        String length = "";
        int count;
        boolean isFlagType;
        switch (type) {
            case UNBOUNDED:
                length = "VAR";
                break;
            case A:
                length = "A";
                break;
            case R:
                length = "R";
                break;
            case G:
                length = "G";
                break;
            case INTEGER: {
                if (headerLine instanceof VCFFormatHeaderLine) {
                    VCFFormatHeaderLine formatHeaderLine = (VCFFormatHeaderLine) headerLine;
                    count = formatHeaderLine.getCount();
                    isFlagType = formatHeaderLine.getType().equals(VCFHeaderLineType.Flag);
                } else {
                    VCFInfoHeaderLine infoHeaderLine = (VCFInfoHeaderLine) headerLine;
                    count = infoHeaderLine.getCount();
                    isFlagType = infoHeaderLine.getType().equals(VCFHeaderLineType.Flag);
                }
                //Weird Flag fields - Number=0 in the VCF header :(
                if (count == 0 && isFlagType) length = "1";
                else length = String.valueOf(count);
                break;
            }
        }
        return length;
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
        throws ParseException {
        String existingVidJson = GenomicsDBUtils.readEntireFile(vidFile);
        GenomicsDBVidMapProto.VidMappingPB.Builder vidMapBuilder = 
                GenomicsDBVidMapProto.VidMappingPB.newBuilder();
        merge(existingVidJson, vidMapBuilder);
        return vidMapBuilder.build();
    }

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

    public static DuplicateVcfFieldNamesCheckResult checkForDuplicateVcfFieldNames(final String vcfFieldName,
	    final int vcfHeaderLineTypeIdx,
	    final String vcfType, final String vcfLengthDescriptor,
	    List<GenomicsDBVidMapProto.GenomicsDBFieldInfo> infoFields,
	    Map<String, List<DuplicateFieldInfoTrackClass>> fieldNameToDuplicateCheckInfo
	    ) {
	List<DuplicateFieldInfoTrackClass> currList = null;
	if(fieldNameToDuplicateCheckInfo.containsKey(vcfFieldName))
	    currList = fieldNameToDuplicateCheckInfo.get(vcfFieldName);
	else {
	    //Why 3? VCF header can have FILTER/INFO/FORMAT fields of the same name
	    List<DuplicateFieldInfoTrackClass> init = new ArrayList(3);
	    for(int i=0;i<3;++i)
		init.add(new DuplicateFieldInfoTrackClass());
	    fieldNameToDuplicateCheckInfo.put(vcfFieldName, init);
	    currList = init;
	}
	assert (currList.size() == 3);
	DuplicateFieldInfoTrackClass currObj = currList.get(vcfHeaderLineTypeIdx);
	if(currObj.isInitialized) {
	    System.err.println("GenomicsDB: "+VidMapExtensions.class.getName()+
		" VCF header has 2 lines with the same field name "+vcfFieldName
		    +" and header line type <"+vcfHeaderLineTypeString[vcfHeaderLineTypeIdx]+">"
		    +" -- only the first line is considered by GenomicsDB");
	    return new DuplicateVcfFieldNamesCheckResult(true, null);
	}
	currObj.vcfType = vcfType;
	currObj.vcfLengthDescriptor = vcfLengthDescriptor;
	//Check if there are other lines of a different type (FILTER, INFO, FORMAT) with the same name
	for(int i=0;i<3;++i) {
	    DuplicateFieldInfoTrackClass otherObj = currList.get(i);
	    if(otherObj.isInitialized) {
		//Check if the other line is of same type and length
		if(otherObj.vcfType.equals(currObj.vcfType) &&
			otherObj.vcfLengthDescriptor.equals(currObj.vcfLengthDescriptor)) {
		    //Simply update the PB object corresponding to the field
		    GenomicsDBVidMapProto.GenomicsDBFieldInfo fieldInfo =
			infoFields.get(otherObj.idxInGenomicsDBFieldInfoList);
		    infoFields.set(otherObj.idxInGenomicsDBFieldInfoList,
			    fieldInfo.toBuilder().
			    addVcfFieldClass(vcfHeaderLineTypeString[vcfHeaderLineTypeIdx])
			    .build());
		    currObj.isInitialized = true;
		    currObj.idxInGenomicsDBFieldInfoList = otherObj.idxInGenomicsDBFieldInfoList;
		    return new DuplicateVcfFieldNamesCheckResult(true, null);
		}
		else { //this is a different type of field, rename
		    currObj.isInitialized = true;
		    currObj.idxInGenomicsDBFieldInfoList = infoFields.size();
		    return new DuplicateVcfFieldNamesCheckResult(false, vcfFieldName+"_"
			    +vcfHeaderLineTypeString[vcfHeaderLineTypeIdx]);
		}
	    }
	}
	currObj.isInitialized = true;
	currObj.idxInGenomicsDBFieldInfoList = infoFields.size();
	//Fields of the same name doesn't exist
	return new DuplicateVcfFieldNamesCheckResult(false, vcfFieldName);
    }
}
