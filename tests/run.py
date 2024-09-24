#!/usr/bin/env python3

#The MIT License (MIT)
#Copyright (c) 2016-2017 Intel Corporation
#Copyright (c) 2019-2020 Omics Data Automation, Inc.
#Copyright (c) 2024 dātma, inc™

#Permission is hereby granted, free of charge, to any person obtaining a copy of 
#this software and associated documentation files (the "Software"), to deal in 
#the Software without restriction, including without limitation the rights to 
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
#the Software, and to permit persons to whom the Software is furnished to do so, 
#subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all 
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
#COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
#IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import json
import tempfile
import subprocess
import hashlib
import os
import sys
import shutil
from collections import OrderedDict
import jsondiff
import errno
import uuid
import tarfile

import common

query_json_template_string="""
{   
        "workspace" : "",
        "array_name" : "",
        "vcf_header_filename" : "inputs/template_vcf_header.vcf",
        "query_column_ranges": [{
            "range_list": [{
                "low": 0,
                "high": 10000000000
            }]
        }],
        "query_row_ranges": [{
            "range_list": [{
                "low": 0,
                "high": 3
            }]
        }],
        "query_filter" : "",
        "reference_genome" : "inputs/chr1_10MB.fasta.gz",
        "attributes" : [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT" ]
}"""

vcf_attributes_order = [ "END", "REF", "ALT", "BaseQRankSum", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "MQ", "RAW_MQ", "MQ0", "DP", "GT", "GQ", "SB", "AD", "PL", "PGT", "PID", "MIN_DP", "DP_FORMAT", "FILTER" ];
asa_vcf_attributes = [ "END", "REF", "ALT", "BaseQRankSum", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "MQ", "RAW_MQ", "MQ0", "DP", "GT", "GQ", "SB", "AD", "PL", "PGT", "PID", "MIN_DP", "DP_FORMAT", "FILTER", "AS_RAW_MQ", "AS_RAW_MQRankSum" ];
attributes_with_DS_ID = [ "REF", "ALT", "BaseQRankSum", "MQ", "RAW_MQ", "MQ0", "ClippingRankSum", "MQRankSum", "ReadPosRankSum", "DP", "GT", "GQ", "SB", "AD", "PL", "DP_FORMAT", "MIN_DP", "PID", "PGT", "DS", "ID" ];
attributes_with_PL_only = [ "PL" ]
attributes_with_MLEAC_only = [ "MLEAC" ]
attributes_with_bigint = [ "big_field1", "big_field2", "AS_big_field3", "big_field4", "AS_big_field5", "DP", "GT" ]
attributes_multi_contig = ["DP", "AF", "GT"]
attributes_multi_contig_empty_contigs = [ "REF", "ALT", "ID", "AD", "DP", "GQ", "GT", "PL", "BaseQRankSum", "HaplotypeScore", "MLEAC", "MLEAF", "MQ", "MQRankSum", "NS", "ReadPosRankSum", "SB", "DP_FORMAT", "AN"] 
# multi contig interval below is missing scaffold_3123 to test that we can deal with non-contiguous contigs when coalescing
multi_contig_empty_contigs_intervals = ['scaffold_3370:1-2985', 'scaffold_3100:1-2123', 'scaffold_3160:1-1884', 'scaffold_3300:1-2080', 'scaffold_3060:1-2290', 'scaffold_3360:1-1235', 'scaffold_3190:1-1726', 'scaffold_3020:1-2530', 'scaffold_3050:1-2317', 'scaffold_3340:1-2888', 'scaffold_2930:1-2674', 'scaffold_3230:1-1594', 'scaffold_3320:1-2574', 'scaffold_3210:1-1863', 'scaffold_3440:1-1034', 'scaffold_3550:1-1938', 'scaffold_3250:1-1485', 'scaffold_3391:1-1147', 'scaffold_3431:1-1055', 'scaffold_3171:1-1840', 'scaffold_3191:1-2934', 'scaffold_3061:1-2903', 'scaffold_2911:1-2734', 'scaffold_2951:1-2618', 'scaffold_3231:1-1591', 'scaffold_2981:1-2544', 'scaffold_3021:1-2419', 'scaffold_3381:1-1181', 'scaffold_3011:1-2971', 'scaffold_3681:1-1022', 'scaffold_2851:1-2950', 'scaffold_3121:1-2043', 'scaffold_3251:1-2376', 'scaffold_2922:1-2707', 'scaffold_3452:1-1858', 'scaffold_3362:1-1482', 'scaffold_3342:1-2686', 'scaffold_3042:1-2357', 'scaffold_3372:1-1208', 'scaffold_3472:1-1331', 'scaffold_3202:1-2247', 'scaffold_3312:1-1324', 'scaffold_3432:1-1048', 'scaffold_3272:1-1426', 'scaffold_3422:1-2139', 'scaffold_3412:1-2813', 'scaffold_3222:1-1606', 'scaffold_3262:1-1448', 'scaffold_3352:1-2879', 'scaffold_3442:1-1032', 'scaffold_3332:1-1283', 'scaffold_3292:1-1613', 'scaffold_3252:1-1483', 'scaffold_3383:1-1162', 'scaffold_3423:1-1189', 'scaffold_3353:1-1243', 'scaffold_2933:1-2669', 'scaffold_2853:1-2935', 'scaffold_2913:1-2832', 'scaffold_2893:1-2797', 'scaffold_3413:1-1105', 'scaffold_3453:1-1003', 'scaffold_3523:1-1130', 'scaffold_3213:1-1653', 'scaffold_3223:1-1705', 'scaffold_3263:1-1970', 'scaffold_3053:1-2603', 'scaffold_3363:1-1230', 'scaffold_3403:1-2418', 'scaffold_3393:1-1141', 'scaffold_3294:1-1369', 'scaffold_3384:1-2287', 'scaffold_3364:1-2412', 'scaffold_3044:1-2336', 'scaffold_3134:1-2599', 'scaffold_3264:1-1441', 'scaffold_2844:1-2966', 'scaffold_3214:1-1641', 'scaffold_3064:1-2391', 'scaffold_3304:1-1992', 'scaffold_3014:1-2633', 'scaffold_3024:1-2408', 'scaffold_3354:1-2537', 'scaffold_3194:1-2061', 'scaffold_2924:1-2704', 'scaffold_3054:1-2300', 'scaffold_3074:1-2965', 'scaffold_3374:1-1200', 'scaffold_3445:1-1028', 'scaffold_3335:1-1924', 'scaffold_3005:1-2485', 'scaffold_2885:1-2936', 'scaffold_2855:1-2931', 'scaffold_3285:1-2448', 'scaffold_3255:1-2564', 'scaffold_3095:1-2137', 'scaffold_3245:1-2263', 'scaffold_3425:1-1075', 'scaffold_2965:1-2576', 'scaffold_3275:1-1417', 'scaffold_3325:1-1298', 'scaffold_3305:1-1337', 'scaffold_3346:1-1260', 'scaffold_3426:1-1074', 'scaffold_3416:1-1094', 'scaffold_3046:1-2329', 'scaffold_3246:1-1501', 'scaffold_3326:1-2466', 'scaffold_3436:1-2815', 'scaffold_3356:1-2921', 'scaffold_3206:1-1669', 'scaffold_3316:1-1314', 'scaffold_3296:1-1362', 'scaffold_3116:1-2637', 'scaffold_3136:1-2710', 'scaffold_3306:1-1470', 'scaffold_3236:1-1544', 'scaffold_3446:1-1027', 'scaffold_3126:1-2028', 'scaffold_3226:1-1920', 'scaffold_3076:1-2214', 'scaffold_3056:1-2299', 'scaffold_3427:1-2622', 'scaffold_2907:1-2749', 'scaffold_3157:1-1898', 'scaffold_3387:1-1156', 'scaffold_2987:1-2534', 'scaffold_3437:1-1040', 'scaffold_3007:1-2593', 'scaffold_2857:1-2929', 'scaffold_3537:1-1011', 'scaffold_3547:1-1702', 'scaffold_3147:1-2230', 'scaffold_3257:1-1463', 'scaffold_3357:1-1239', 'scaffold_3107:1-2869', 'scaffold_3388:1-1155', 'scaffold_3368:1-1220', 'scaffold_3068:1-2254', 'scaffold_2838:1-2988', 'scaffold_3008:1-2908', 'scaffold_3058:1-2824', 'scaffold_3648:1-2796', 'scaffold_3328:1-1293', 'scaffold_3398:1-1131', 'scaffold_3528:1-2635', 'scaffold_3348:1-1258', 'scaffold_2968:1-2573', 'scaffold_3478:1-2690', 'scaffold_3228:1-1597', 'scaffold_2868:1-2983', 'scaffold_3088:1-2159', 'scaffold_3298:1-2443', 'scaffold_3448:1-1020', 'scaffold_3239:1-1965', 'scaffold_2879:1-2846', 'scaffold_2919:1-2717', 'scaffold_3269:1-1431', 'scaffold_3129:1-2114', 'scaffold_3389:1-1153', 'scaffold_3179:1-2709', 'scaffold_3379:1-1639', 'scaffold_3139:1-1964', 'scaffold_3259:1-1665', 'scaffold_3289:1-2781', 'scaffold_2889:1-2826']
default_segment_size = 40

def create_query_json(ws_dir, test_name, query_param_dict):
    test_dict=json.loads(query_json_template_string);
    test_dict["workspace"] = ws_dir
    test_dict["array_name"] = test_name
    if('query_column_ranges' in query_param_dict):
        test_dict["query_column_ranges"] = query_param_dict["query_column_ranges"]
    else:
        test_dict['scan_full'] = True
    if("vid_mapping_file" in query_param_dict):
        test_dict["vid_mapping_file"] = query_param_dict["vid_mapping_file"];
    if("callset_mapping_file" in query_param_dict):
        test_dict["callset_mapping_file"] = query_param_dict["callset_mapping_file"];
    if("attributes" in query_param_dict):
        test_dict["attributes"] = query_param_dict["attributes"]
    if('segment_size' in query_param_dict):
        test_dict['segment_size'] = query_param_dict['segment_size'];
    else:
        test_dict['segment_size'] = default_segment_size;
    if('produce_GT_field' in query_param_dict):
        test_dict['produce_GT_field'] = query_param_dict['produce_GT_field'];
    if('produce_FILTER_field' in query_param_dict):
        test_dict['produce_FILTER_field'] = query_param_dict['produce_FILTER_field'];
    if('sites_only_query' in query_param_dict):
        test_dict['sites_only_query'] = query_param_dict['sites_only_query']
    if('produce_GT_with_min_PL_value_for_spanning_deletions' in query_param_dict):
        test_dict['produce_GT_with_min_PL_value_for_spanning_deletions'] = \
                query_param_dict['produce_GT_with_min_PL_value_for_spanning_deletions']
    if('query_sample_names' in query_param_dict):
        test_dict['query_sample_names'] = query_param_dict['query_sample_names']
        if('query_row_ranges' in test_dict):
            del test_dict['query_row_ranges']
    if('query_filter' in query_param_dict):
        test_dict['query_filter'] = query_param_dict['query_filter']
    if('max_genotype_count' in query_param_dict):
        test_dict['max_genotype_count'] = query_param_dict['max_genotype_count']
    if('reference_genome' in query_param_dict):
        test_dict['reference_genome'] =  query_param_dict['reference_genome']
    return test_dict;


loader_json_template_string="""
{
    "row_based_partitioning" : false,
    "column_partitions" : [
        {"begin": 0, "workspace":"", "array_name": "" }
    ],
    "callset_mapping_file" : "",
    "vid_mapping_file" : "inputs/vid.json",
    "size_per_column_partition": 700 ,
    "treat_deletions_as_intervals" : true,
    "vcf_header_filename": "inputs/template_vcf_header.vcf",
    "reference_genome" : "inputs/chr1_10MB.fasta.gz",
    "num_parallel_vcf_files" : 1,
    "do_ping_pong_buffering" : false,
    "offload_vcf_output_processing" : false,
    "discard_vcf_index": true,
    "produce_combined_vcf": true,
    "produce_tiledb_array" : true,
    "delete_and_create_tiledb_array" : true,
    "compress_tiledb_array" : false,
    "segment_size" : 1048576,
    "num_cells_per_tile" : 3
}""";

def create_loader_json(ws_dir, test_name, test_params_dict):
    test_dict=json.loads(loader_json_template_string);
    if('column_partitions' in test_params_dict):
        test_dict['column_partitions'] = test_params_dict['column_partitions'];
    test_dict["column_partitions"][0]["workspace"] = ws_dir;
    test_dict["column_partitions"][0]["array_name"] = test_name;
    test_dict["callset_mapping_file"] = test_params_dict['callset_mapping_file'];
    if('vid_mapping_file' in test_params_dict):
        test_dict['vid_mapping_file'] = test_params_dict['vid_mapping_file'];
    if('size_per_column_partition' in test_params_dict):
        test_dict['size_per_column_partition'] = test_params_dict['size_per_column_partition'];
    if('segment_size' in test_params_dict):
        test_dict['segment_size'] = test_params_dict['segment_size'];
    else:
        test_dict['segment_size'] = default_segment_size;
    if('reference_genome' in test_params_dict):
        test_dict['reference_genome'] =  test_params_dict['reference_genome']
    return test_dict;

def run_cmd(cmd, expected_to_pass, errstring='Sanity check : ', tmpdir=None, print_if_error=True,
        dont_capture_output=False):
    print(f"cmd={cmd}")
    if(dont_capture_output):
        returncode = subprocess.call(cmd, shell=True)
        stdout_string = None
        stderr_string = None
    else:
        pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_string, stderr_string = pid.communicate()
        returncode = pid.returncode
    is_success = (returncode == 0)
    if(not is_success):
        if(print_if_error):
            sys.stderr.write(errstring+ ' -- failed command -- '+cmd+'\n')
            if(not dont_capture_output):
                sys.stderr.write('STDOUT:\n')
                sys.stderr.write(stdout_string.decode('utf-8'))
                sys.stderr.write('STDERR:\n')
                sys.stderr.write(stderr_string.decode('utf-8'))
        if(expected_to_pass):
            cleanup_and_exit(None, -1)
    return (is_success, stdout_string, stderr_string)

def get_file_content_and_md5sum(filename):
    with open(filename, 'rb') as fptr:
        data = fptr.read();
        md5sum_hash_str = str(hashlib.md5(data).hexdigest())
        fptr.close();
        return (data, md5sum_hash_str);

def print_diff(golden_output, test_output):
    print("=======Golden output:=======");
    print(golden_output);
    print("=======Test output:=======");
    print(test_output);
    print("=======END=======");

def modify_query_column_ranges_for_PB(test_query_dict):
    if('query_column_ranges' in test_query_dict):
        original_query_column_ranges = test_query_dict['query_column_ranges']
        new_query_column_ranges = []
        for curr_entry in original_query_column_ranges:
            if(type(curr_entry) is dict and 'range_list' in curr_entry):
                new_interval_list = []
                for curr_interval in curr_entry['range_list']:
                    if(type(curr_interval) is dict and 'low' in curr_interval
                            and 'high' in curr_interval):
                        new_interval_list.append({'column_interval': { 'tiledb_column_interval':
                            { 'begin': curr_interval['low'], 'end': curr_interval['high'] } } })
                new_entry = { 'column_or_interval_list': new_interval_list }
                new_query_column_ranges.append(new_entry)
        if not new_query_column_ranges:
            test_query_dict['query_column_ranges'] = original_query_column_ranges
        else:
            test_query_dict['query_column_ranges'] = new_query_column_ranges

def bcftools_compare(bcftools_path, exe_path, outfilename, outfilename_golden, outfile_string, golden_string):
    with open(outfilename, "w") as out_fd:
        out_fd.write(outfile_string.decode('utf-8'))
    bcf_cmd = bcftools_path+' view -Oz -o '+outfilename+'.gz '+outfilename+' && '+bcftools_path+' index -f -t '+outfilename+'.gz'
    success,bcf_stdout,stderr_string = run_cmd(bcf_cmd, False, 'Call to bcftools failed ')
    if(not success):
        return True

    with open(outfilename_golden, "w") as out_golden_fd:
        out_golden_fd.write(golden_string)
    bcf_cmd = bcftools_path+' view -Oz -o '+outfilename_golden+'.gz '+outfilename_golden+' && '+bcftools_path+' index -f -t '+outfilename_golden+'.gz'
    success,bcf_stdout,stderr_string = run_cmd(bcf_cmd, False, 'Call to bcftools failed ')
    if(not success):
        return True

    vcfdiff_cmd = exe_path+os.path.sep+'vcfdiff '+outfilename+'.gz '+outfilename_golden+'.gz'
    success, vcfdiff_stdout,stderr_string = run_cmd(vcfdiff_cmd, False, 'vcfdiff cmd failed ')
    if(not success):
        return True
    
    if(len(vcfdiff_stdout) == 0):
        return False
    else:
        sys.stderr.write('vcfdiff check failed. Cmd: '+vcfdiff_cmd+'\n')
        sys.stderr.write('Returned: '+vcfdiff_stdout+'\n')
        return True

def create_json_file(tmpdir, test_name, query_type, test_dict):
    query_json_filename = tmpdir+os.path.sep+test_name+ (('_'+query_type) if query_type else '') +'.json'
    with open(query_json_filename, 'w') as fptr:
        json.dump(test_dict, fptr, indent=4, separators=(',', ': '))
        fptr.close()
    return query_json_filename

def cleanup_and_exit(tmpdir, exit_code):
    if(exit_code == 0 and tmpdir):
       shutil.rmtree(tmpdir, ignore_errors=True)
    sys.exit(exit_code);

def substitute_workspace_dir(jsonfile, wsdir):
    import fileinput
    for line in fileinput.input(jsonfile, inplace=True):
        print(line.replace("#WORKSPACE_DIR#", wsdir)),

def test_with_java_options(tmpdir, lib_path, jacoco):
    sys.stdout.write("Testing with Java Options...\n")
    compatDir = tmpdir+'/compat_100_test'
    tar = tarfile.open('inputs/compatibility.pre.1.0.0.test.tar')
    tar.extractall(path=compatDir)
    loader_json = compatDir+'/t0_1_2.json'
    query_json = compatDir+'/t0_1_2_java_vcf.json'
    substitute_workspace_dir(loader_json, compatDir)
    substitute_workspace_dir(query_json, compatDir)

    # Sanity Test java options supported by GenomicsDB
    query_command = 'java'+jacoco+' -Dgenomicsdb.library.path='+lib_path+' -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -ea TestGenomicsDB --query -l '+loader_json+' '+query_json
    query_java_options_success, stdout_string, stderr_string = run_cmd(query_command, True,
        'Query with java options failed ')
    actual_md5sum = str(hashlib.md5(stdout_string).hexdigest())
    golden_str, golden_md5sum = get_file_content_and_md5sum('golden_outputs/java_t0_1_2_vcf_at_0')
    if (actual_md5sum != golden_md5sum):
        sys.stderr.write('Query with java options : '+query_command+' failed. Results did not match with golden output.\n');
        cleanup_and_exit(tmpdir, -1)

    # Test if exception stack trace is logged when there is no tiledb_workspace for example
    os.remove(compatDir+'/ws/__tiledb_workspace.tdb')
    query_java_options_success, stdout_string, stderr_string = run_cmd(query_command, False,
        '', print_if_error=False)
    if (str(stderr_string).find('Native Stack Trace:') < 0):
        sys.stderr.write('Query with java options did not throw expected exception: '+query_command+' failed\n')
        cleanup_and_exit(None, -1)
    sys.stdout.write("Successful\n")

def test_pre_1_0_0_query_compatibility(tmpdir):
    sys.stdout.write("Testing compatibility with workspace/arrays from releases before 1.0.0...") 
    compatDir = tmpdir+'/compat_100_test'
    tar = tarfile.open('inputs/compatibility.pre.1.0.0.test.tar')
    tar.extractall(path=compatDir)
    loader_json = compatDir+'/t0_1_2.json'
    query_json = compatDir+'/t0_1_2_java_vcf.json'
    substitute_workspace_dir(loader_json, compatDir)
    substitute_workspace_dir(query_json, compatDir)
    query_command = 'java -ea TestGenomicsDB --query -l '+loader_json+' '+query_json
    success, stdout_string, stderr_string = run_cmd(query_command, True, 'Compatibility Test Query failed ')
    actual_md5sum = str(hashlib.md5(stdout_string).hexdigest())
    golden_str, golden_md5sum = get_file_content_and_md5sum('golden_outputs/java_t0_1_2_vcf_at_0')
    if (actual_md5sum != golden_md5sum):
        sys.stderr.write('Compatibility Test Query : '+query_command+' failed. Results did not match with golden output.\n');
        cleanup_and_exit(tmpdir, -1)
    sys.stdout.write("Successful\n")

def test_query_compatibility_with_old_schema_verion_1(tmpdir):
    sys.stdout.write("Testing query compatibility with arrays created with old array schema version 1...") 
    compatDir = tmpdir+'/compat_schema_version_1_test'
    tar = tarfile.open('inputs/compatibility.array.schema.version.1.tar')
    tar.extractall(path=compatDir)
    loader_json = compatDir+'/t0_1_2.json'
    query_json = compatDir+'/t0_1_2_query.json'
    substitute_workspace_dir(loader_json, compatDir)
    substitute_workspace_dir(query_json, compatDir)
    query_command = 'java -ea TestGenomicsDB --query -l '+loader_json+' '+query_json
    success, stdout_string,stderr_string = run_cmd(query_command, True, 'Compatibility Test Query failed ')
    actual_md5sum = str(hashlib.md5(stdout_string).hexdigest())
    golden_str, golden_md5sum = get_file_content_and_md5sum('golden_outputs/java_t0_1_2_vcf_sites_only_at_0')
    if (actual_md5sum != golden_md5sum):
        sys.stderr.write('Compatibility Test Query : '+query_command+' failed. Results did not match with golden output.\n');
        cleanup_and_exit(tmpdir, -1)
    sys.stdout.write("Successful\n")

def tool_sanity_checks(exe_path):
    gt_mpi_gather_path = exe_path+os.path.sep+'gt_mpi_gather';
    try:
        st = os.stat(exe_path)
    except os.error:
        sys.stderr.write("Could not find " + gt_mpi_gather_path + '\n');
        sys.exit(-1);
    run_cmd(gt_mpi_gather_path, False, print_if_error=False);
    run_cmd(gt_mpi_gather_path + ' --help', True);
    run_cmd(gt_mpi_gather_path + ' --h', True);
    run_cmd(gt_mpi_gather_path + ' --version', True);
    run_cmd(gt_mpi_gather_path + ' --gibberish', False, print_if_error=False);
    run_cmd(gt_mpi_gather_path + ' -j non_existent_query_json', False, print_if_error=False);
    run_cmd(gt_mpi_gather_path + ' -j inputs/tools/non_existent_callset.json', False, print_if_error=False);
    run_cmd(gt_mpi_gather_path + ' -j inputs/tools/non_existent_vidmap.json', False, print_if_error=False);
    run_cmd(gt_mpi_gather_path + ' -j inputs/tools/non_existent_workspace_array.json', False, print_if_error=False);
    run_cmd(gt_mpi_gather_path + ' -j inputs/tools/unparseable_query.json', False, print_if_error=False);

def test_pb_configs(tmpdir, ctest_dir):
    pb_query_test_configs = [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                        }],
                      "query_row_ranges": [{
                          "range_list": [{
                              "low": 0,
                              "high": 1
                              },
                              {
                              "low": 3,
                              "high": 3
                              }
                              ]
                          }],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'produce_GT_field': True,
                      'segment_size': 100,
                    },
                    { "query_samples_names": [
                        "HG00141", "NA12878"
                        ],
                       "query_contig_intervals": [
                           { "contig": "1" },
                           { "contig": "1", "begin": 100 },
                           { "contig": "1", "begin": 100, "end": 500 },
                           ],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'sites_only_query': True,
                      'segment_size': 100,
                    },
                ]
    for test_pb_dict in pb_query_test_configs:
        query_dict = create_query_json('/tmp/ws', 'pb_test', test_pb_dict)
        if('scan_full' in query_dict):
            del query_dict['scan_full']
        modify_query_column_ranges_for_PB(query_dict)
        query_json_filename = create_json_file(tmpdir, 'pb_test', None, query_dict)
        cmd = ctest_dir+'/ctests pb_query_config_test --query-json-file '+query_json_filename
        run_cmd(cmd, True, 'PB query configs test')

def main():
    #Initial Setup
    if (len(sys.argv) < 3):
        sys.stderr.write('Usage ./run.py <build_dir> <install_dir> [<build_type>] [<split_batch_num>]\n')
        sys.stderr.write('   Optional Argument 3 - build_type=Release|Coverage|...\n')
        sys.stderr.write('   Optional Argument 4 - Tests are batched into two to help with Travis builds on MacOS\n')
        sys.stderr.write('                         specify batch with split_batch_num=1|2\n')
        sys.stderr.write('                         Otherwise all tests are run.\n\n')
        sys.exit(-1);
    build_dir=os.path.abspath(sys.argv[1])
    install_dir=os.path.abspath(sys.argv[2])
    ctest_dir = os.path.join(build_dir, 'src/test/cpp')
    exe_path = os.path.join(install_dir,'bin')
    lib_path = os.path.join(install_dir,'lib')
    tool_sanity_checks(exe_path);
    if (len(sys.argv) >= 4):
        build_type = sys.argv[3]
    else:
        build_type = "default"
    #Switch to tests directory
    parent_dir=os.path.dirname(os.path.realpath(__file__))
    os.chdir(parent_dir)
    top_tmpdir = tempfile.mkdtemp()
    bcftools_path = shutil.which("bcftools")
    java_import_dir = top_tmpdir+os.path.sep + 'java_import'
    os.mkdir(java_import_dir)
    test_pb_configs(top_tmpdir, ctest_dir)
    common.setup_classpath(build_dir)
    jacoco, jacoco_report_cmd = common.setup_jacoco(build_dir, build_type)
    loader_tests0 = [
            { "name" : "t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                    "query_sample_names" : [ "HG00141", "HG01958" ] ,
                    "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_vcf_at_0",
                        } },
                    { "query_column_ranges": [{ 'column_or_interval_list': [{'column': {'tiledb_column': 12000}}, 
                                                                            {'column': {'tiledb_column': 12142}},
                                                                            {'column': {'tiledb_column': 12144}},
                                                                            {'column': {'tiledb_column': 12160}},
                                                                            {'column': {'tiledb_column': 12290}},
                                                                            {'column': {'tiledb_column': 12294}},
                                                                            {'column': {'tiledb_column': 14000}},
                                                                            {'column': {'tiledb_column': 17384}},
                                                                            {'column': {'tiledb_column': 18000}}]}],
                       "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_multiple_positions",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_multiple_positions",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_multiple_positions",
                        } },
                    { "query_column_ranges": [{ 'column_or_interval_list': [{'column_interval': {'tiledb_column_interval': {'begin': 0, 'end': 1000000}}}]}], 
                      "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges" : [{
                        "range_list": [{
                             "low": 0,
                             "high": 1000000000
                          }]
                      }],
                      "sites_only_query": True,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_sites_only_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_sites_only_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        },{
                            "low": 12141,
                            "high": 12141
                        }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100_12141",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        },{
                            "low": 12141,
                            "high": 12141
                        },{
                            "low": 12150,
                            "high": 12150
                        }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100_12141_12150",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        },
                            {
                                "low": 12141,
                                "high": 12150
                            }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100_12141_to_12150",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12100,
                            "high": 12100
                        },
                            {
                                "low": 12141,
                                "high": 12150
                            },
                            {
                                "low": 12300,
                                "high": 12300
                            },
                            {
                                "low": 17384,
                                "high": 17384
                            }]
                    }], "golden_output": {   #
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12100_12141_to_12150_12300_17384",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "attributes": attributes_with_PL_only,
                        "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0_with_PL_only",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],#vid and callset jsons passed through query json
                        "query_without_loader": True,
                        "vid_mapping_file": "inputs/vid.json",
                        "callset_mapping_file": "inputs/callsets/t0_1_2.json",
                        "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "produce_FILTER_field": True, "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0_with_FILTER",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig",
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'chromosome_intervals': [ '1:1-12160', '1:12161-12200', '1:12201-18000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                'generate_array_name_from_partition_bounds': True,
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 18000
                        }]
                    }],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig_vcf_0_18000",
                        } },
                    {
                        "query_contig_interval": {
                            "contig": "1",
                            "begin": 12151,
                            "end": 18000
                        },
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig_vcf_12150_18000",
                        } }
                ]
            },
             { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig_query_non_existent_contig",
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'chromosome_intervals': [ '1:1-12160', '1:12161-12200', '1:12201-18000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                'generate_array_name_from_partition_bounds': True,
                "query_params": [
                    {
                        "query_contig_interval": {
                            "contig": "2",
                            "begin": 12151,
                            "end": 18000
                        },
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_empty_vcf",
                        } }
                ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_incremental_t0_1_2_multi_contig",
                'callset_mapping_file': 'inputs/callsets/t0_1_2.0.json',
                'callset_mapping_file1': 'inputs/callsets/t0_1_2.1.json',
                'chromosome_intervals': [ '1:1-12160', '1:12161-12200', '1:12201-18000' ],
                'consolidate_only': '3',
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                'generate_array_name_from_partition_bounds': True,
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 18000
                        }]
                    }],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig_vcf_0_18000",
                        "consolidate_java_vcf"   : "golden_outputs/java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig_vcf_0_18000",
                        } },
                    {
                        "query_contig_interval": {
                            "contig": "1",
                            "begin": 12151,
                            "end": 18000
                        },
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_genomicsdb_importer_from_vcfs_t0_1_2_multi_contig_vcf_12150_18000",
                        } }
                ]
            },
            { "name" : "t0_1_2_filter", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_csv.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                    'query_filter' : 'DP > 100',
                    "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0_with_DP_filter",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0_with_DP_filter",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0_with_DP_filter",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0_with_DP_filter",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0_with_DP_filter",
                        } },
                    ]
            },
            { "name" : "t0_1_2_csv", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_csv.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "t0_overlapping", 'golden_output': 'golden_outputs/t0_overlapping',
                'callset_mapping_file': 'inputs/callsets/t0_overlapping.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12202,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "vcf"        : "golden_outputs/t0_overlapping_at_12202",
                        }
                    }
                ]
            },
            { "name" : "t0_overlapping_at_12202", 'golden_output': 'golden_outputs/t0_overlapping_at_12202',
                'callset_mapping_file': 'inputs/callsets/t0_overlapping.json',
                'column_partitions': [ {"begin": 12202, "workspace":"", "array_name": "" }]
            },
            { "name" : "spanning_with_partition_straddle", 'golden_output': 'golden_outputs/spanning_with_partition_straddle_load_stdout',
                'callset_mapping_file': 'inputs/callsets/spanning_partition.json',
                'column_partitions': [ {"begin": 1274365, "workspace":"", "array_name": ""}]
            },
            { "name" : "deletion_with_partition_straddle", 'golden_output': 'golden_outputs/deletion_with_partition_straddle_load_stdout',
                'callset_mapping_file': 'inputs/callsets/t6.json',
                'column_partitions': [ {"begin": 8029500, "workspace":"", "array_name": ""}]
            },
    ];
    loader_tests1 = [
            { "name" : "t6_7_8", 'golden_output' : 'golden_outputs/t6_7_8_loading',
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_0",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_0",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_0",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                       "max_genotype_count": 8,
                       "golden_output": {
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_0_max_genotype_count_8",
                        } },
                    { "query_column_ranges" : [{
                        "range_list": [{
                          "low": 0,
                          "high": 1000000000
                        }]
                       }],
                      "sites_only_query": True, "golden_output": {
                        "vcf"        : "golden_outputs/t6_7_8_vcf_sites_only_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 8029500,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_8029500",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_8029500",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_8029500",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_8029500",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 8029500,
                            "high": 8029500
                        }]
                    }], "golden_output": {
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_8029500-8029500",
                        } }
                    ]
            },
            { "name" : "t6_SNV_7_MNV_8", 'golden_output' : 'golden_outputs/t6_SNV_t7_MNV_t8_loading.vcf',
                'callset_mapping_file': 'inputs/callsets/t6_SNV_7_MNV_8.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "vcf"        : "golden_outputs/t6_SNV_t7_MNV_t8_loading.vcf",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t6_7_8_multi_contig",
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                'chromosome_intervals': [ '1:1-8029500','1:8029501-8029501', '1:8029502-10000000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                'generate_array_name_from_partition_bounds': True,
                "query_params": [
                    {   'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_0",
                        } },
                    {
                        "query_contig_interval": {
                            "contig": "1",
                            "begin": 8029501,
                            "end": 8029510
                        },
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_8029500",
                        } },
                    {
                        "query_contig_interval": {
                            "contig": "1",
                            "begin": 8029502,
                            "end": 8029502
                        },
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_8029501",
                        } }
                ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_incremental_t6_7_8_multi_contig",
                'callset_mapping_file': 'inputs/callsets/t6_7_8.0.json',
                'callset_mapping_file1': 'inputs/callsets/t6_7_8.1.json',
                'chromosome_intervals': [ '1:1-8029500','1:8029501-8029501', '1:8029502-10000000' ],
                'consolidate_only': '2',
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                'generate_array_name_from_partition_bounds': True,
                "query_params": [
                    {   'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_0",
                        "consolidate_java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_0",
                        } },
                    {
                        "query_contig_interval": {
                            "contig": "1",
                            "begin": 8029501,
                            "end": 8029510
                        },
                        'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_8029500",
                        } },
                    {
                        "query_contig_interval": {
                            "contig": "1",
                            "begin": 8029502,
                            "end": 8029502
                        },
                        'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "java_vcf"   : "golden_outputs/java_t6_7_8_vcf_at_8029501",
                        } }
                ]
            },
            { "name" : "java_t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0_phased_GT",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150_phased_GT",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150_phased_GT",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "java_buffer_stream_t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_buffer.json',
                'stream_name_to_filename_mapping': 'inputs/callsets/t0_1_2_buffer_mapping.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "java_buffer_stream_multi_contig_t0_1_2", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_buffer.json',
                'stream_name_to_filename_mapping': 'inputs/callsets/t0_1_2_buffer_mapping.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_0",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_calls_at_12150",
                        "variants"   : "golden_outputs/t0_1_2_variants_at_12150",
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "test_new_fields", 'golden_output' : 'golden_outputs/t6_7_8_new_field_gatk.vcf',
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                'vid_mapping_file': 'inputs/vid_MLEAC_MLEAF.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      "attributes" : attributes_with_MLEAC_only, "golden_output": {
                        "calls"        : "golden_outputs/test_new_fields_MLEAC_only.json",
                        } },
                    ]
            },
            { "name" : "test_info_combine_ops0", 'golden_output' : 'golden_outputs/info_ops0.vcf',
                'callset_mapping_file': 'inputs/callsets/info_ops.json',
                'vid_mapping_file': 'inputs/vid_info_ops0.json'
            },
            { "name" : "test_info_combine_ops1", 'golden_output' : 'golden_outputs/info_ops1.vcf',
                'callset_mapping_file': 'inputs/callsets/info_ops.json',
                'vid_mapping_file': 'inputs/vid_info_ops1.json'
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2",
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_incremental_t0_1_2",
                'callset_mapping_file': 'inputs/callsets/t0_1_2.0.json',
                'callset_mapping_file1': 'inputs/callsets/t0_1_2.1.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_0",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_0",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 12150,
                            "high": 1000000000
                        }]
                    }],
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_vcf_at_12150",
                        "batched_vcf": "golden_outputs/t0_1_2_vcf_at_12150",
                        "java_vcf"   : "golden_outputs/java_t0_1_2_vcf_at_12150",
                        } }
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t6_7_8",
                'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_0_phased_GT",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_0",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 8029500,
                            "high": 1000000000
                        }]
                    }],
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_8029500_phased_GT",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_8029500_phased_GT",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_8029500",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_8029500",
                        } }
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_incremental_t6_7_8",
                'callset_mapping_file': 'inputs/callsets/t6_7_8.0.json',
                'callset_mapping_file1': 'inputs/callsets/t6_7_8.1.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_0_phased_GT",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_0",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_0",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 8029500,
                            "high": 1000000000
                        }]
                    }],
                        "vid_mapping_file": "inputs/vid_phased_GT.json",
                        'callset_mapping_file': 'inputs/callsets/t6_7_8.json',
                        "golden_output": {
                        "calls"      : "golden_outputs/t6_7_8_calls_at_8029500_phased_GT",
                        "variants"   : "golden_outputs/t6_7_8_variants_at_8029500_phased_GT",
                        "vcf"        : "golden_outputs/t6_7_8_vcf_at_8029500",
                        "batched_vcf": "golden_outputs/t6_7_8_vcf_at_8029500",
                        } }
                    ]
            },
            { "name" : "t0_1_2_combined", 'golden_output' : 'golden_outputs/t0_1_2_combined',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_combined.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "vcf"        : "golden_outputs/t0_1_2_combined",
                        "batched_vcf": "golden_outputs/t0_1_2_combined",
                        } },
                    ]
            }, 
            { "name" : "test_flag_field", 'golden_output' : 'golden_outputs/t0_1_2_DS_ID_vcf_at_0',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'vid_mapping_file': 'inputs/vid_DS_ID.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "attributes": attributes_with_DS_ID, "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_DS_ID_calls_at_0",
                        "variants"   : "golden_outputs/t0_1_2_DS_ID_variants_at_0",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2_with_DS_ID",
                'vid_mapping_file': 'inputs/vid_DS_ID_phased_GT.json',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        'vid_mapping_file': 'inputs/vid_DS_ID_phased_GT.json',
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "attributes": attributes_with_DS_ID, "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_DS_ID_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t0_1_2_DS_ID_variants_at_0_phased_GT",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_incremental_t0_1_2_with_DS_ID",
                'vid_mapping_file': 'inputs/vid_DS_ID_phased_GT.json',
                'callset_mapping_file': 'inputs/callsets/t0_1_2.0.json',
                'callset_mapping_file1': 'inputs/callsets/t0_1_2.1.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        'vid_mapping_file': 'inputs/vid_DS_ID_phased_GT.json',
                        'callset_mapping_file': 'inputs/callsets/t0_1_2.json',
                        "attributes": attributes_with_DS_ID, "golden_output": {
                        "calls"      : "golden_outputs/t0_1_2_DS_ID_calls_at_0_phased_GT",
                        "variants"   : "golden_outputs/t0_1_2_DS_ID_variants_at_0_phased_GT",
                        } },
                    ]
            },
            { "name" : "t0_1_2_as_array", 'golden_output' : 'golden_outputs/t0_1_2_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_as_array.json',
                "vid_mapping_file": "inputs/vid_as_array.json",
            },
            { "name" : "t0_with_missing_PL_SB_fields", 'golden_output' : 'golden_outputs/t0_with_missing_PL_SB_fields_t1.vcf',
                'callset_mapping_file': 'inputs/callsets/t0_with_missing_PL_SB_fields_t1.json',
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }], "golden_output": {
                        "calls"      : "golden_outputs/t0_with_missing_PL_SB_fields_t1_calls.json",
                        } },
                    ]
            },
            { "name" : "t0_haploid_triploid_1_2_3_triploid_deletion",
                'golden_output' : 'golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_loading',
                'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                'size_per_column_partition': 1200,
                'segment_size': 100,
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'segment_size': 100,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                        }],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'produce_GT_field': True,
                      'segment_size': 100,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf_produce_GT",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf_produce_GT",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                        }],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'produce_GT_field': True,
                      'produce_GT_with_min_PL_value_for_spanning_deletions': True,
                      'segment_size': 100,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf_produce_GT_for_min_value_PL",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf_produce_GT_for_min_PL",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                        }],
                      'callset_mapping_file': 'inputs/callsets/t0_haploid_triploid_1_2_3_triploid_deletion.json',
                      "vid_mapping_file": "inputs/vid_DS_ID_phased_GT.json",
                      'sites_only_query': True,
                      'segment_size': 100,
                      "golden_output": {
                        "vcf"        : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_vcf_sites_only",
                        "java_vcf"   : "golden_outputs/t0_haploid_triploid_1_2_3_triploid_deletion_java_vcf_sites_only",
                        } },
                ]
            },
            { "name" : "t0_1_2_all_asa", 'golden_output' : 'golden_outputs/t0_1_2_all_asa_loading',
                'callset_mapping_file': 'inputs/callsets/t0_1_2_all_asa.json',
                'vid_mapping_file': 'inputs/vid_all_asa.json',
                'size_per_column_partition': 3000,
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      "force_override": True,
                      'segment_size': 100,
                      "attributes": asa_vcf_attributes,
                        "golden_output": {
                        "vcf"      : "golden_outputs/t0_1_2_all_asa_loading",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2_all_asa",
                'callset_mapping_file': 'inputs/callsets/t0_1_2_all_asa.json',
                'vid_mapping_file': 'inputs/vid_all_asa.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "force_override": True,
                        'segment_size': 100,
                        "attributes": asa_vcf_attributes,
                        "golden_output": {
                        "vcf"      : "golden_outputs/t0_1_2_all_asa_loading",
                        "java_vcf"   : "golden_outputs/t0_1_2_all_asa_java_query_vcf",
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_t0_1_2_no_remap_missing",
                'callset_mapping_file': 'inputs/callsets/t0_1_2_all_asa.json',
                'vid_mapping_file': 'inputs/vid_no_remap_missing.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "force_override": True,
                        'segment_size': 100,
                        "attributes": asa_vcf_attributes,
                        "golden_output": {
                        "vcf"      : "golden_outputs/t0_1_2_loading_no_remap_missing"
			# no java_vcf here because htsjdk has a bug that causes issues when AD/PL have missing values 
			# https://github.com/broadinstitute/gatk/issues/6744#issuecomment-674975007
                        } },
                    ]
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_incremental_t0_1_2_all_asa",
                'callset_mapping_file': 'inputs/callsets/t0_1_2_all_asa.0.json',
                'callset_mapping_file1': 'inputs/callsets/t0_1_2_all_asa.1.json',
                'vid_mapping_file': 'inputs/vid_all_asa.json',
                'chromosome_intervals': [ '1:1-100000000' ],
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                        "force_override": True,
                        'segment_size': 100,
                        "attributes": asa_vcf_attributes,
                        "golden_output": {
                        "vcf"      : "golden_outputs/t0_1_2_all_asa_loading",
                        "java_vcf"   : "golden_outputs/t0_1_2_all_asa_java_query_vcf",
                        } },
                    ]
            },
            { "name" : "min_PL_spanning_deletion", 'golden_output' : 'golden_outputs/min_PL_spanning_deletion_load_stdout',
                'callset_mapping_file': 'inputs/callsets/min_PL_spanning_deletion.json',
                "vid_mapping_file": "inputs/vid_phased_GT.json",
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      'produce_GT_field': True, "golden_output": {
                        "vcf"        : "golden_outputs/min_PL_spanning_deletion_vcf_no_min_PL",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      'produce_GT_field': True,
                      'produce_GT_with_min_PL_value_for_spanning_deletions': True,
                      "golden_output": {
                        "vcf"        : "golden_outputs/min_PL_spanning_deletion_vcf",
                        } }
                    ]
            },
            { "name" : "t6_7_8_asa", 'golden_output' : 'golden_outputs/t6_7_8_asa.vcf',
                'callset_mapping_file': 'inputs/callsets/t6_7_8_asa.json',
                'vid_mapping_file': 'inputs/vid_all_asa.json',
                'size_per_column_partition': 3000,
                "query_params": [
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 0,
                            "high": 1000000000
                        }]
                    }],
                      "force_override": True,
                      'segment_size': 100,
                      "attributes": asa_vcf_attributes,
                        "golden_output": {
                        "vcf"      : "golden_outputs/t6_7_8_asa.vcf",
                        } },
                    { "query_column_ranges": [{
                        "range_list": [{
                            "low": 8029500,
                            "high": 1000000000
                        }]
                    }],
                      "force_override": True,
                      'segment_size': 100,
                      "attributes": asa_vcf_attributes,
                        "golden_output": {
                        "vcf"      : "golden_outputs/t6_7_8_asa_vcf_at_8029500",
                        } },
                    ]
            },
            { "name" : "test_info_combine_bigint", 'golden_output' : 'golden_outputs/bigint.vcf',
              'callset_mapping_file': 'inputs/callsets/bigint.json',
              'vid_mapping_file': 'inputs/vid_bigint.json',
              'size_per_column_partition': 4096,
              "query_params": [
                  {
                    'attributes': attributes_with_bigint,
                    "force_override": True,
                    'segment_size': 4096,
                    'pass_as_vcf': True,
                    "query_column_ranges": [{
                      "range_list": [{
                          "low": 0,
                          "high": 1000000000
                      }]
                  }], "golden_output": {
                      "java_vcf"   : "golden_outputs/bigint_java.vcf",
                      } },
               ],
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_multi_contig_coalesce",
              'callset_mapping_file': 'inputs/callsets/multi_contig.json',
              'vid_mapping_file': 'inputs/vid_multi_contig.json',
              'reference_genome': 'inputs/chr1_2_3_4_70Kb.fasta.gz',
              'size_per_column_partition': 4096,
              'generate_array_name_from_partition_bounds': True,
              'coalesce_to_num_partitions': '2',
              'chromosome_intervals': [ '1:1-249250621', '2:1-243199373', '3:1-198022430', '4:1-191154276' ],
              "query_params": [
                  {
                    'attributes': attributes_multi_contig,
                    "force_override": True,
                    'segment_size': 4096,
                    'pass_as_vcf': True,
                    'produce_GT_field': True,
                    'reference_genome': 'inputs/chr1_2_3_4_70Kb.fasta.gz',
                    "query_column_ranges": [{
                      "range_list": [{
                          "low": 0,
                          "high": 1000000000
                      }]
                  }], "golden_output": {
                      "java_vcf"   : "golden_outputs/multi_contig_java_vcf",
                      } },
               ],
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_incremental_multi_contig_coalesce",
              'callset_mapping_file': 'inputs/callsets/multi_contig.0.json',
              'callset_mapping_file1': 'inputs/callsets/multi_contig.1.json',
              'vid_mapping_file': 'inputs/vid_multi_contig.json',
              'reference_genome': 'inputs/chr1_2_3_4_70Kb.fasta.gz',
              'size_per_column_partition': 4096,
              'generate_array_name_from_partition_bounds': True,
              'coalesce_to_num_partitions': '2',
              'consolidate_only': '0',
              'chromosome_intervals': [ '1:1-249250621', '2:1-243199373', '3:1-198022430', '4:1-191154276' ],
              "query_params": [
                  {
                    'attributes': attributes_multi_contig,
                    "force_override": True,
                    'segment_size': 4096,
                    'pass_as_vcf': True,
                    'produce_GT_field': True,
                    'reference_genome': 'inputs/chr1_2_3_4_70Kb.fasta.gz',
                    "query_column_ranges": [{
                      "range_list": [{
                          "low": 0,
                          "high": 1000000000
                      }]
                  }], "golden_output": {
                      "java_vcf"   : "golden_outputs/multi_contig_java_vcf",
                      "consolidate_java_vcf"   : "golden_outputs/multi_contig_java_vcf",
                      } },
               ],
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_multi_contig_empty_contigs_coalesce",
              'callset_mapping_file': 'inputs/callsets/multi_contig_empty_contigs.json',
              'vid_mapping_file': 'inputs/vid_multi_contig_empty_contigs.json',
              'reference_genome': 'inputs/Pop_tri_nonchrom_subset.fa.gz',
              'size_per_column_partition': 4096,
              'generate_array_name_from_partition_bounds': True,
              'coalesce_to_num_partitions': '7',
              'consolidate_only': '1',
              'chromosome_intervals': multi_contig_empty_contigs_intervals,
              "query_params": [
                  {
                    'attributes': attributes_multi_contig_empty_contigs,
                    "force_override": True,
                    'segment_size': 4096,
                    'pass_as_vcf': True,
                    'produce_GT_field': True,
                    'reference_genome': 'inputs/Pop_tri_nonchrom_subset.fa.gz',
                    "query_column_ranges": [{
                      "range_list": [{
                          "low": 0,
                          "high": 1000000000
                      }]
                  }], "golden_output": {
                      "java_vcf"   : "golden_outputs/multi_contig_empty_contigs_java_vcf",
                      "consolidate_java_vcf"   : "golden_outputs/multi_contig_empty_contigs_java_vcf",
                      } },
               ],
            },
            { "name" : "java_genomicsdb_importer_from_vcfs_duplicate_fields",
              'callset_mapping_file': 'inputs/callsets/duplicate_fields.json',
              'vid_mapping_file': 'inputs/vid_multi_contig_empty_contigs.json',
              'size_per_column_partition': 4096,
              'generate_array_name_from_partition_bounds': True,
              'chromosome_intervals': [ '1:1-150121160' ],
              "query_params": [
                  {
                    'attributes': attributes_multi_contig,
                    "force_override": True,
                    'segment_size': 4096,
                    'pass_as_vcf': True,
                    'produce_GT_field': True,
                    "query_column_ranges": [{
                      "range_list": [{
                          "low": 0,
                          "high": 1000000000
                      }]
                  }], "golden_output": {
                      "java_vcf"   : "golden_outputs/duplicate_fields_java_vcf",
                      } },
               ],
            },
    ];
    if(len(sys.argv) < 5):
        loader_tests = loader_tests0 + loader_tests1
    elif(sys.argv[4] == '1'):
        loader_tests = loader_tests0
    else:
        loader_tests = loader_tests1
    test_names_set = set()
    for test_params_dict in loader_tests:
        test_name = test_params_dict['name']
        if(test_name in test_names_set):
          new_name = test_name+'_'+str(uuid.uuid4())
          sys.stderr.write('Duplicate test names in loader_tests '+test_name+' renamed to '+new_name+'\n')
          test_name = new_name
        test_names_set.add(test_name)
        #New tmpdir per test
        tmpdir = top_tmpdir + os.path.sep + test_name
        try:
          os.makedirs(tmpdir, 0o755)
        except os.error as e:
          if(e.errno == errno.EEXIST):
            pass
          else:
            raise e
        ws_dir = tmpdir + os.path.sep + 'ws'
        if(test_name.find('java_genomicsdb_importer_from_vcfs') != -1):
            ws_dir = java_import_dir+os.path.sep+test_name
        if not os.path.isdir(ws_dir):
            create_workspace_cmd = exe_path+os.path.sep+'create_genomicsdb_workspace '+ws_dir
            run_cmd(create_workspace_cmd, True)
        if(test_name.find('java_genomicsdb_importer_from_vcfs_duplicate_fields') != -1):
            test_params_dict['vid_mapping_file'] = tmpdir+os.path.sep+'vid.json'
        test_loader_dict = create_loader_json(ws_dir, test_name, test_params_dict);
        incremental_load = False
        if(test_name == "t0_1_2"):
            test_loader_dict["compress_tiledb_array"] = True;
        # loader json won't get used with genomicsdb importer, but will need to be accurate 
        # for queries (specifically --produce-Broad-GVCF). Point to the combined callset for queries
        if(test_name.find('java_genomicsdb_importer_from_vcfs_incremental') != -1):
            incr_callset = test_loader_dict["callset_mapping_file"]
            test_loader_dict["callset_mapping_file"] = incr_callset.replace('.0', '')
        loader_json_filename = create_json_file(tmpdir, 'loader', None, test_loader_dict)
        if(test_name  == 'java_t0_1_2'):
            import_cmd = 'java'+jacoco+' -ea TestGenomicsDB --load '+loader_json_filename

        elif(test_name == 'java_buffer_stream_multi_contig_t0_1_2'):
            import_cmd = 'java'+jacoco+' -ea TestBufferStreamGenomicsDBImporter -iterators '+loader_json_filename+' ' + \
                    test_params_dict['stream_name_to_filename_mapping'] + \
                    ' 1024 0 0 100 true '
        elif(test_name == 'java_buffer_stream_t0_1_2'):
            import_cmd = 'java'+jacoco+' -ea TestBufferStreamGenomicsDBImporter '+loader_json_filename \
                    +' '+test_params_dict['stream_name_to_filename_mapping']
        elif(test_name.find('java_genomicsdb_importer_from_vcfs') != -1):
            arg_list = ''
            for interval in test_params_dict['chromosome_intervals']:
                arg_list += ' -L '+interval
            arg_list += ' -w ' + ws_dir +' --use_samples_in_order ' + ' --batchsize=1 '
            arg_list += ' --vidmap-output '+ tmpdir + os.path.sep + 'vid.json'
            arg_list += ' --callset-output '+ tmpdir + os.path.sep + 'callsets.json'
            if('generate_array_name_from_partition_bounds' not in test_params_dict or
                    not test_params_dict['generate_array_name_from_partition_bounds']):
                arg_list += ' -A ' + test_name
            file_list = ''
            count=0
            with open(test_params_dict['callset_mapping_file'], 'rb') as cs_fptr:
                callset_mapping_dict = json.load(cs_fptr, object_pairs_hook=OrderedDict)
                for callset_name, callset_info in callset_mapping_dict['callsets'].items():
                    file_list += ' '+callset_info['filename'];
                    count=count+1
                cs_fptr.close();

            coalesce_arg = ''
            if('coalesce_to_num_partitions' in test_params_dict):
                coalesce_arg = ' --coalesce-multiple-contigs '+test_params_dict['coalesce_to_num_partitions']
            import_cmd = 'java'+jacoco+' -ea TestGenomicsDBImporterWithMergedVCFHeader --size_per_column_partition 16384 ' \
                                   '--segment_size 10485760'+arg_list+file_list
            if('consolidate_only' in test_params_dict):
                import_cmd_consolidate = import_cmd + ' --consolidate-only '+test_params_dict['consolidate_only']
            import_cmd = import_cmd + coalesce_arg
            if(test_name.find('incremental') != -1):
                incremental_load = True
                file_list = ''
                with open(test_params_dict['callset_mapping_file1'], 'rb') as cs_fptr:
                    callset_mapping_dict = json.load(cs_fptr, object_pairs_hook=OrderedDict)
                    for callset_name, callset_info in callset_mapping_dict['callsets'].items():
                        file_list += ' '+callset_info['filename'];
                cs_fptr.close();
                import_cmd_incremental = 'java'+jacoco+' -ea TestGenomicsDBImporterWithMergedVCFHeader --size_per_column_partition 16384 ' \
                                               '--segment_size 10485760 --incremental '+str(count)+arg_list+file_list+coalesce_arg
        else:
            import_cmd = exe_path+os.path.sep+'vcf2genomicsdb '+loader_json_filename
        
        loader_cmd_success,stdout_string,stderr_string = run_cmd(import_cmd, True, 'Loader test: '+test_name)
        md5sum_hash_str = str(hashlib.md5(stdout_string).hexdigest())
        if('golden_output' in test_params_dict):
            golden_stdout, golden_md5sum = get_file_content_and_md5sum(test_params_dict['golden_output']);
            if(golden_md5sum != md5sum_hash_str):
                is_error = True
                if(bcftools_path is None):
                    sys.stderr.write('Did not find bcftools in path. If you are able to add that to the path, we can get a more precise diff\n')
                else:
                    outfilename = tmpdir+os.path.sep+test_name+'_loader.vcf'
                    outfilename_golden = tmpdir+os.path.sep+test_name+'_loader_golden'+'.vcf'
                    is_error = (bcftools_compare(bcftools_path, exe_path, outfilename, outfilename_golden, stdout_string, golden_stdout))
                if(is_error):
                    sys.stderr.write('Loader stdout mismatch for test: '+test_name+'\n');
                    print_diff(golden_stdout, stdout_string);
                    cleanup_and_exit(None, -1)

        if incremental_load:
            inc_load_success,stdout_string,stderr_string = run_cmd(import_cmd_incremental, True, 'Loader test: '+test_name)

        if('query_params' in test_params_dict):
            for query_param_dict in test_params_dict['query_params']:
                test_query_dict = create_query_json(ws_dir, test_name, query_param_dict)
                if(test_name.find('java_genomicsdb_importer_from_vcfs') != -1 and
                        'generate_array_name_from_partition_bounds' in test_params_dict
                        and test_params_dict['generate_array_name_from_partition_bounds']):
                    if('array' in test_query_dict):
                        del test_query_dict['array']
                    if('array_name' in test_query_dict):
                        del test_query_dict['array_name']
                    if('query_column_ranges' in test_query_dict):
                        del test_query_dict['query_column_ranges']
                        test_query_dict['scan_full'] = True
                query_types_list = [
                        ('calls','--print-calls'),
                        ('variants',''),
                        ('vcf','--produce-Broad-GVCF -p 128'),
                        ('batched_vcf','--produce-Broad-GVCF -p 128'),
                        ('java_vcf', ''),
                        ('consolidate_java_vcf', ''),
                        ('consolidate_and_vcf', '--produce-Broad-GVCF'), #keep as the last query test
                        ]
                for query_type,cmd_line_param in query_types_list:
                    if('golden_output' in query_param_dict and query_type in query_param_dict['golden_output']):
                        if((query_type == 'vcf' or query_type == 'batched_vcf' or query_type.find('java_vcf') != -1)
                                and 'force_override' not in query_param_dict):
                            test_query_dict['attributes'] = vcf_attributes_order;
                        if(query_type.find('java_vcf') != -1):
                            modify_query_column_ranges_for_PB(test_query_dict)
                        query_json_filename = create_json_file(tmpdir, 'q', query_type, test_query_dict)
                        if(query_type.find('java_vcf') != -1):
                            loader_argument = loader_json_filename;
                            misc_args = ''
                            if("query_without_loader" in query_param_dict and query_param_dict["query_without_loader"]):
                                loader_argument = '""'
                            if('query_contig_interval' in query_param_dict):
                                query_contig_interval_dict = query_param_dict['query_contig_interval']
                                misc_args += ('--chromosome '+query_contig_interval_dict['contig'] \
                                        + ' --begin %d --end %d')%(query_contig_interval_dict['begin'],
                                                query_contig_interval_dict['end'])
                            if('pass_as_vcf' in query_param_dict and query_param_dict['pass_as_vcf']):
                                misc_args += ' --pass_as_vcf '
                            if query_type == "consolidate_java_vcf":
                                run_cmd(import_cmd_consolidate, True, 'Consolidate only: {} '.format(test_name))
                            query_command = 'java'+jacoco+' -ea TestGenomicsDB --query -l '+loader_argument+' '+query_json_filename \
                                + ' ' + misc_args;
                        else:
                            if(query_type == 'consolidate_and_vcf'):
                                run_cmd(exe_path+os.path.sep+'consolidate_genomicsdb_array '+ws_dir+' '+test_name,
                                    True, 'TileDB array consolidation failed for test '+test_name+' : ')
                            loader_argument = ' -l '+loader_json_filename;
                            if("query_without_loader" in query_param_dict and query_param_dict["query_without_loader"]):
                                loader_argument = ''
                            query_command = (exe_path+os.path.sep+'gt_mpi_gather -s %d'+loader_argument
                                + ' -j '
                                +query_json_filename+' '+cmd_line_param)%(test_query_dict['segment_size']);
                        query_success,stdout_string,stderr_string = run_cmd(query_command, True,
                                'Query test: {}-{} '.format(test_name, query_type))
                        md5sum_hash_str = str(hashlib.md5(stdout_string).hexdigest())
                        golden_stdout, golden_md5sum = get_file_content_and_md5sum(query_param_dict['golden_output'][query_type]);
                        if(golden_md5sum != md5sum_hash_str):
                            is_error = True;
                            print("Command="+query_command+"\n");
                            #do JSON diff for variant and call format print
                            json_diff_result = None
                            if(query_type in set(['calls', 'variants'])):
                                try:
                                    golden_stdout_dict = json.loads(golden_stdout);
                                    test_stdout_dict = json.loads(stdout_string);
                                    json_diff_result = jsondiff.diff(golden_stdout_dict, test_stdout_dict);
                                    if(len(json_diff_result) == 0):
                                        is_error = False;
                                except:
                                    json_diff_result = None;
                                    is_error = True
                            # for VCF format try to use bcftools and call vcfdiff
                            else:
                                if(bcftools_path is None):
                                    sys.stderr.write('Did not find bcftools in path. If you are able to add that to the path, we can get a more precise diff\n')
                                else:
                                    outfilename = tmpdir+os.path.sep+test_name+'_'+query_type+'.vcf'
                                    outfilename_golden = tmpdir+os.path.sep+test_name+'_'+query_type+'_golden'+'.vcf'
                                    is_error = bcftools_compare(bcftools_path, exe_path, outfilename, outfilename_golden, stdout_string, golden_stdout)
                            if(is_error):
                                sys.stderr.write('Mismatch in query test: '+test_name+'-'+query_type+'\n');
                                sys.stderr.write('Query command: '+query_command+'\n')
                                print_diff(golden_stdout, stdout_string);
                                if(json_diff_result):
                                    print(json.dumps(json_diff_result, indent=4, separators=(',', ': ')));
                                cleanup_and_exit(None, -1);
                        #ignore some corner cases for now
                        if(query_type == 'vcf'
                            and 'query_filter' not in query_param_dict
                            ):
                            #These tests fail because PL remapping is WIP
                            skip_GT_matching = (test_name.find('t0_haploid_triploid_1_2_3_triploid_deletion') != -1
                              or test_name.find('min_PL_spanning_deletion') != -1)
                            cmd = ctest_dir+'/ctests --durations yes columnar_gvcf_iterator_test ' \
                                    + ' --query-json-file ' + query_json_filename \
                                    + ' --golden-output-file '+query_param_dict['golden_output'][query_type] \
                                    + ' --loader-json-file '+loader_json_filename \
                                    + (' --skip-GT-matching' if(skip_GT_matching) else '')
                            run_cmd(cmd, True, 'Ctests '+test_name, dont_capture_output=True)
    test_with_java_options(top_tmpdir, lib_path, jacoco)
    test_pre_1_0_0_query_compatibility(top_tmpdir)
    test_query_compatibility_with_old_schema_verion_1(top_tmpdir)
    rc = common.report_jacoco_coverage(jacoco_report_cmd)
    cleanup_and_exit(top_tmpdir, rc)

if __name__ == '__main__':
    main()
