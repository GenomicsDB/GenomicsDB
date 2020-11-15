/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

import org.genomicsdb.exception.GenomicsDBException;
import org.genomicsdb.importer.GenomicsDBImporter;
import org.genomicsdb.model.CommandLineImportConfig;
import org.json.simple.parser.ParseException;

import java.util.Arrays;
import java.io.IOException;

public final class TestGenomicsDBImporterWithMergedVCFHeader {

  public static void main(String[] args) throws IOException, GenomicsDBException, ParseException, InterruptedException {
    // this expects --coalesce-multiple-contigs to be the last argument, if it exists
    // similarly --consolidate-only. and they should be mututally exclusive
    int numPartitions = 0, doConsolidate = -1;
    if (args[args.length-2].equals("--coalesce-multiple-contigs")) {
      numPartitions = Integer.parseInt(args[args.length-1]);
      args = Arrays.copyOf(args, args.length-2);
    }
    else if (args[args.length-2].equals("--consolidate-only")) {
      doConsolidate = Integer.parseInt(args[args.length-1]);
      args = Arrays.copyOf(args, args.length-2);
    }
    CommandLineImportConfig config = new CommandLineImportConfig("TestGenomicsDBImporterWithMergedVCFHeader", args);
    GenomicsDBImporter importer = new GenomicsDBImporter(config);
    if (numPartitions > 0) {
      importer.coalesceContigsIntoNumPartitions(numPartitions);
    }
    if (doConsolidate >= 0) {
      importer.doConsolidate(doConsolidate);
    }
    else{
      importer.executeImport();
    }
  }

}
