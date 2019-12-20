/**
 * The MIT License (MIT) Copyright (c) 2019 Omics Data Automation
 *
 * <p>Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 * and associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * <p>The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * <p>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 * BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
import gnu.getopt.Getopt;
import gnu.getopt.LongOpt;
import org.apache.spark.SparkContext;
import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.DataFrameReader;
import org.apache.spark.sql.SparkSession;
import org.apache.spark.sql.types.*;
import org.apache.spark.sql.types.DataTypes.*;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import org.genomicsdb.model.GenomicsDBExportConfiguration;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.nio.charset.Charset;
import java.util.*;

import static org.apache.spark.sql.functions.bround;
import static org.apache.spark.sql.functions.col;

import com.googlecode.protobuf.format.JsonFormat;

public final class TestGenomicsDBDataSourceV2 {

  public static String[] floatFields = {
    "BaseQRankSum", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "MQ", "RAW_MQ", "MLEAF"
  };

  private static String readFile(String path, Charset encoding) throws IOException {
    byte[] encoded = Files.readAllBytes(Paths.get(path));
    return new String(encoded, encoding);
  }

  public static boolean schemaContainsField(StructType schema, String s) {
    String[] fields = schema.fieldNames();
    for (String fs : fields) {
      if (fs.equals(s)) {
        return true;
      }
    }
    return false;
  }

  public static StructType getSchemaFromQuery(File query, String vidmapping)
      throws IOException, ParseException {
    try {
      JSONParser qparser = new JSONParser();
      FileReader queryJsonReader = new FileReader(query);
      JSONObject obj = null;
      JSONParser vparser = new JSONParser();
      FileReader vidJsonReader = new FileReader(vidmapping);
      JSONObject vobj = null;
      try {
        obj = (JSONObject) qparser.parse(queryJsonReader);
        vobj = (JSONObject) vparser.parse(vidJsonReader);
      } catch (ParseException | IOException e) {
        queryJsonReader.close();
        vidJsonReader.close();
        throw e;
      }

      /* we open both vid and query fields
       * and iterate through all attributes in query json
       * matching them up with vid fields, and creating a schema
       */
      JSONObject fields = (JSONObject) vobj.get("fields");
      JSONArray attrArray = (JSONArray) obj.get("query_attributes");
      ArrayList<StructField> schemaArray = new ArrayList<>();
      if (attrArray != null) {
        for (Object attrObj : attrArray) {
          String attr = (String) attrObj;
          if (attr.equals("DP_FORMAT")) {
            schemaArray.add(new StructField(attr, DataTypes.NullType, true, Metadata.empty()));
          } else if (!(attr.equals("GT")
              || attr.equals("REF")
              || attr.equals("FILTER")
              || attr.equals("ALT")
              || attr.equals("ID"))) {
            JSONObject v = (JSONObject) fields.get(attr);
            JSONArray fieldClass = (JSONArray) v.get("vcf_field_class");
            // these will all be nulltype, and the datasourcev2 api for genomicsdb
            // takes care of correctly assigning the datatype from VID mapping
            if (fieldClass != null) {
              schemaArray.add(new StructField(attr, DataTypes.NullType, true, Metadata.empty()));
            }
          }
        }
      }

      queryJsonReader.close();
      StructType st = new StructType(schemaArray.toArray(new StructField[0]));
      return st;
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      return null;
    }
  }

  public static void main(final String[] args) throws IOException,
         org.json.simple.parser.ParseException {
    LongOpt[] longopts = new LongOpt[6];
    longopts[0] = new LongOpt("loader", LongOpt.REQUIRED_ARGUMENT, null, 'l');
    longopts[1] = new LongOpt("query", LongOpt.REQUIRED_ARGUMENT, null, 'q');
    longopts[2] = new LongOpt("vid", LongOpt.REQUIRED_ARGUMENT, null, 'v');
    longopts[3] = new LongOpt("hostfile", LongOpt.REQUIRED_ARGUMENT, null, 'h');
    longopts[4] = new LongOpt("spark_master", LongOpt.REQUIRED_ARGUMENT, null, 's');
    longopts[5] = new LongOpt("use-query-protobuf", LongOpt.NO_ARGUMENT, null, 'p');

    if (args.length < 8) {
      System.err.println(
          "Usage:\n\t--loader <loader.json> --query <query.json> --vid <vid.json> "
              + "--spark_master <sparkMaster>"
              +"\nOptional args:\n--hostfile <hostfile> --use-query-protobuf");
      System.exit(-1);
    }
    String loaderFile, queryFile, hostfile, vidMapping, sparkMaster, jarDir;
    boolean useQueryProtobuf = false;
    loaderFile = queryFile = hostfile = sparkMaster = vidMapping = "";
    Getopt g = new Getopt("TestGenomicsDBSparkHDFS", args, "l:q:h:s:v:p", longopts);
    int c = -1;
    String optarg;

    while ((c = g.getopt()) != -1) {
      switch (c) {
        case 'l':
          loaderFile = g.getOptarg();
          break;
        case 'q':
          queryFile = g.getOptarg();
          break;
        case 'v':
          vidMapping = g.getOptarg();
          break;
        case 'h':
          hostfile = g.getOptarg();
          break;
        case 's':
          sparkMaster = g.getOptarg();
          break;
        case 'p':
          useQueryProtobuf = true;
          break;
        default:
          System.err.println("Unknown command line option " + g.getOptarg());
          System.exit(-1);
      }
    }
    SparkSession spark = SparkSession.builder().appName("TestGenomicsDBDataSourceV2").getOrCreate();

    Path dstdir = Paths.get("").toAbsolutePath();
    Path qSrc = Paths.get(queryFile);
    Path lSrc = Paths.get(loaderFile);
    File qDstFile = null;
    if (!useQueryProtobuf) {
      qDstFile = File.createTempFile("query", ".json", new File(dstdir.toString()));
      qDstFile.deleteOnExit();
    }
    File lDstFile = File.createTempFile("loader", ".json", new File(dstdir.toString()));
    lDstFile.deleteOnExit();
    if (!useQueryProtobuf) {
      Files.copy(qSrc, qDstFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
    }
    Files.copy(lSrc, lDstFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

    SparkContext sc = spark.sparkContext();
    if (!useQueryProtobuf) {
      sc.addFile(qDstFile.getName());
    }
    sc.addFile(lDstFile.getName());

    StructType schema = null;
    try {
      schema = getSchemaFromQuery(new File(queryFile), vidMapping);
    } catch (ParseException e) {
      e.printStackTrace();
    }

    Dataset<Row> variants;
    DataFrameReader reader = spark.read()
            .format("org.genomicsdb.spark.GenomicsDBDataSourceV2")
            .schema(schema)
            .option("genomicsdb.input.loaderjsonfile", lDstFile.getName());
    if (!hostfile.isEmpty()) {
      reader = reader.option("genomicsdb.input.mpi.hostfile", hostfile);
    }
    if (useQueryProtobuf) {
      GenomicsDBExportConfiguration.ExportConfiguration.Builder builder =
      GenomicsDBExportConfiguration.ExportConfiguration.newBuilder();
      String jsonString = readFile(queryFile, Charset.defaultCharset());
      JsonFormat.merge(jsonString, builder);
      byte[] pb = builder.build().toByteArray();
      String pbString = Base64.getEncoder().encodeToString(pb);

      reader = reader.option("genomicsdb.input.queryprotobuf", pbString);
    } else {
      reader = reader.option("genomicsdb.input.queryjsonfile", qDstFile.getName());
    }
    variants = reader.load();

    String tempDir = "./" + UUID.randomUUID().toString();

    // change number format for our floats so they're easier to read/compare
    for (String s : floatFields) {
      if (schemaContainsField(schema, s)) {
        variants = variants.withColumn(s, bround(col(s), 3));
      }
    }

    variants.coalesce(1).orderBy("contig", "startPos").drop("variantType").write().json(tempDir);

    File f = new File(tempDir);
    File[] matchingFiles =
        f.listFiles(
            new FilenameFilter() {
              public boolean accept(File dir, String name) {
                return name.endsWith("json");
              }
            });

    BufferedReader br = new BufferedReader(new FileReader(matchingFiles[0]));
    String line;
    StringJoiner sj = new StringJoiner(",", "[", "]");
    while ((line = br.readLine()) != null) {
      sj.add(line);
    }
    System.out.println(sj.toString());
    br.close();

    Files.walk(Paths.get(tempDir))
        .map(Path::toFile)
        .sorted((o1, o2) -> -o1.compareTo(o2))
        .forEach(File::delete);
  }
}
