/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2018 Intel Corporation
 * Copyright (c) 2018-2019 Omics Data Automation, Inc.
 * Copyright (c) 2024 dātma, inc™
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

package org.genomicsdb.reader;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.variant.bcf2.BCF2Codec;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.genomicsdb.importer.extensions.JsonFileExtensions;
import org.genomicsdb.model.GenomicsDBExportConfiguration;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static org.genomicsdb.Constants.CHROMOSOME_FOLDER_DELIMITER_SYMBOL_REGEX;
import static org.genomicsdb.GenomicsDBUtils.getArrayColumnBounds;

/**
 * Iterator over {@link htsjdk.variant.variantcontext.VariantContext} objects.
 * Uses {@link GenomicsDBQueryStream} to obtain combined gVCF records
 * (as BCF2) from TileDB/GenomicsDB
 */
public class GenomicsDBFeatureIterator<T extends Feature, SOURCE> implements CloseableTribbleIterator<T>, JsonFileExtensions {

    private static Logger logger = LogManager.getLogger(GenomicsDBFeatureIterator.class);

    private class GenomicsDBQueryStreamParamsHolder {

        public String loaderJSONFile;
        public GenomicsDBExportConfiguration.ExportConfiguration queryPB;
        public String contig;
        public long begin;
        public long end;

        GenomicsDBQueryStreamParamsHolder(final String loaderJSONFile, 
                final GenomicsDBExportConfiguration.ExportConfiguration queryPB,
                final String contig, final long begin, final long end) {
            this.loaderJSONFile = loaderJSONFile;
            this.queryPB = queryPB;
            this.contig = contig;
            this.begin = begin;
            this.end = end;
        }
    };

    private FeatureCodecHeader featureCodecHeader;
    private FeatureCodec<T, SOURCE> codec;
    private List<GenomicsDBQueryStreamParamsHolder> queryParamsList;
    private GenomicsDBTimer timer;
    private boolean closedBefore;
    private SOURCE currentSource;
    private int currentIndexInQueryParamsList;

    /**
     * Constructor
     *
     * @param loaderJSONFile     GenomicsDB loader JSON configuration file
     * @param queryPB            GenomicsDB query protobuf object
     * @param arrayNames         List of array names
     * @param featureCodecHeader htsjdk Feature codec header
     * @param codec              FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                           and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @throws IOException       when data cannot be read from the stream
     */
    GenomicsDBFeatureIterator(final String loaderJSONFile, 
            final GenomicsDBExportConfiguration.ExportConfiguration queryPB,
            final Optional<List<String>> arrayNames, 
            final FeatureCodecHeader featureCodecHeader, 
            final FeatureCodec<T, SOURCE> codec) throws IOException {
        this(loaderJSONFile, queryPB, arrayNames, featureCodecHeader, codec, "", OptionalInt.empty(),
                OptionalInt.empty());
    }

    /**
     * Constructor
     *
     * @param loaderJSONFile     GenomicsDB loader JSON configuration file
     * @param queryPB            GenomicsDB query protobuf
     * @param arrayNames         List of array names
     * @param featureCodecHeader htsjdk Feature codec header
     * @param codec              FeatureCodec, currently only {@link htsjdk.variant.bcf2.BCF2Codec}
     *                           and {@link htsjdk.variant.vcf.VCFCodec} are tested
     * @param chr                contig name
     * @param start              start position (1-based)
     * @param end                end position, inclusive (1-based)
     * @throws IOException       when data cannot be read from the stream
     */
    GenomicsDBFeatureIterator(final String loaderJSONFile, 
            final GenomicsDBExportConfiguration.ExportConfiguration queryPB,
            final Optional<List<String>> arrayNames,
            final FeatureCodecHeader featureCodecHeader, final FeatureCodec<T, SOURCE> codec,
            final String chr, final OptionalInt start, final OptionalInt end) throws IOException {
        this.featureCodecHeader = featureCodecHeader;
        this.codec = codec;
        if (arrayNames.isPresent()) {
            this.queryParamsList = arrayNames.orElse(Collections.emptyList()).stream().map(array -> {
                GenomicsDBExportConfiguration.ExportConfiguration newQueryPB =
                    GenomicsDBExportConfiguration.ExportConfiguration.newBuilder(queryPB)
                    .setArrayName(array).build();
                GenomicsDBQueryStreamParamsHolder params;
                String[] ref = array.split(CHROMOSOME_FOLDER_DELIMITER_SYMBOL_REGEX);
                if(ref.length == 3) {
                    int iStart = Integer.parseInt(ref[1]);
                    int iEnd = Integer.parseInt(ref[2]);
                    if (start.isPresent() && end.isPresent()) {
                        params = new GenomicsDBQueryStreamParamsHolder(loaderJSONFile, 
                                newQueryPB, ref[0], Math.max(start.getAsInt(), iStart), 
                                Math.min(end.getAsInt(), iEnd));
                    }
                    else if (!start.isPresent() && !end.isPresent()) {
                        params = new GenomicsDBQueryStreamParamsHolder(loaderJSONFile, 
                                newQueryPB, ref[0], iStart, iEnd);
                    }
                    else {
                        throw new RuntimeException("Start and End must either be both specified or both unspecified");
                    }
                }
                else {
                    if (start.isPresent() && end.isPresent()) {
                        params = new GenomicsDBQueryStreamParamsHolder(loaderJSONFile,
                                newQueryPB, chr, start.getAsInt(), end.getAsInt());
                    }
                    else {
                        long[] bounds = getArrayColumnBounds(newQueryPB.getWorkspace(), array);
                        params = new GenomicsDBQueryStreamParamsHolder(loaderJSONFile,
                                newQueryPB, "", bounds[0], bounds[1]);
                    }
                }

                return params;
            }).collect(Collectors.toList());
        }
        else {
            boolean isChrEmpty = chr.isEmpty();
            this.queryParamsList = Collections.singletonList(
                    new GenomicsDBQueryStreamParamsHolder(loaderJSONFile, queryPB, chr,
                    start.isPresent() ? start.getAsInt() : isChrEmpty ? 0 : 1,
                    end.isPresent() ? end.getAsInt() : isChrEmpty? 0 : Integer.MAX_VALUE));
        }
        if (queryParamsList.isEmpty()) {
            String message = "GenomicsDB workspace does not have sources for input query interval: " + chr;
            if (start.isPresent()) {
                message += ":" + start.getAsInt();
                if (end.isPresent()) {
                    message += "-" + end.getAsInt();
                }
            }
            logger.warn(message);
        }
        this.currentIndexInQueryParamsList = -1;
        this.currentSource = null;
        setNextSourceAsCurrent();
        this.timer = new GenomicsDBTimer();
        this.closedBefore = false;
    }

  @Override
  public boolean hasNext() {
    // While loop since the next source might not return any data, but subsequent sources might
    while ((this.currentSource == null || this.codec.isDone(this.currentSource))
        && this.currentIndexInQueryParamsList < this.queryParamsList.size()) {
      try {
        setNextSourceAsCurrent();
      } catch (IOException e) {
        throw new RuntimeException(e.getLocalizedMessage());
      }
    }
    boolean isDone = (this.currentSource == null || this.codec.isDone(this.currentSource));
    if (isDone) close();
    return !isDone;
  }

    @Override
    public T next() {
        try {
            this.timer.start();
            if (this.codec.isDone(this.currentSource)) throw new RuntimeException("No more valid records exist");
            T nextObj = this.codec.decode(this.currentSource);
            this.timer.stop();
            return nextObj;
        } catch (IOException e) {
            throw new RuntimeException("Unknown codec exception");
        }
    }

    @Override
    public void close() {
        if (!this.closedBefore) {
            this.timer.print("GenomicsDB iterator next() timer", System.err);
            this.codec.close(this.currentSource);
            this.closedBefore = true;
        }
    }

    @Override
    public Iterator<T> iterator() {
        return this;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported in Iterators");
    }

    private void setNextSourceAsCurrent() throws IOException {
        if(this.currentSource != null)
            this.codec.close(this.currentSource);
        ++(this.currentIndexInQueryParamsList);
        if (this.currentIndexInQueryParamsList < this.queryParamsList.size()) {
            boolean readAsBCF = this.codec instanceof BCF2Codec;
            GenomicsDBQueryStreamParamsHolder currParams = this.queryParamsList.get(this.currentIndexInQueryParamsList);
            GenomicsDBQueryStream queryStream = new GenomicsDBQueryStream(currParams.loaderJSONFile, currParams.queryPB,
                    currParams.contig, currParams.begin, currParams.end, readAsBCF);
            this.currentSource = this.codec.makeSourceFromStream(queryStream);
            if (readAsBCF) { //BCF2 codec provides size of header
                long numByteToSkip = this.featureCodecHeader.getHeaderEnd();
                long numBytesSkipped = queryStream.skip(numByteToSkip);
                if(numBytesSkipped != numByteToSkip)
                    throw new IOException("Could not skip header in GenomicsDBQueryStream - header is "
                            +numByteToSkip+" bytes long but skip() could only bypass "+numBytesSkipped+" bytes");
            }
            else //VCF Codec must parse out header again since getHeaderEnd() returns 0
                this.codec.readHeader(this.currentSource); //no need to store header anywhere
        }
    }
}
