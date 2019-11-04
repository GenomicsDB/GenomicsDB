#!/usr/bin/env python

#The MIT License (MIT)
#Copyright (c) 2019 Omics Data Automation, Inc.

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

import sys
import json
from collections import OrderedDict
from google.protobuf import json_format

#GenomicsDB modules
from genomicsdb.ga4gh import common_pb2
from genomicsdb.ga4gh import variants_pb2

#Ignore the next 5 lines - hack to read in output of gt_mpi_gather --print-GA4GH-calls
top_dict = json.load(sys.stdin, object_pairs_hook=OrderedDict)
for v_dict in top_dict['variants']:
  ga4gh_variant = variants_pb2.Variant()
  json_format.Parse(json.dumps(v_dict), ga4gh_variant)
  print('BEGIN_NEW_CALL')
  #From here, the code walks through the GA4GH structure and prints out data
  #Consider this as an example to traverse the structure
  print('Position %s:%d-%d'%(ga4gh_variant.reference_name, ga4gh_variant.start, ga4gh_variant.end))
  print('REF %s\nALT: %s'%(ga4gh_variant.reference_bases, ','.join(ga4gh_variant.alternate_bases)))
  #currently the following loop will have only 1 element
  for ga4gh_call in ga4gh_variant.calls:
    print('Sample %s'%ga4gh_call.call_set_name)
    print('GT %s'%('/'.join(ga4gh_call.genotype)))
    #other attributes
    for attribute_name, ga4gh_attribute_values_list_wrapper in ga4gh_call.attributes.attr.iteritems():
      ga4gh_attribute_values_list = ga4gh_attribute_values_list_wrapper.values
      if(attribute_name == 'BaseQRankSum'):
	assert len(ga4gh_attribute_values_list) == 1
	ga4gh_attribute_value = ga4gh_attribute_values_list[0]
	assert ga4gh_attribute_value.HasField('double_value')
	print('BaseQRankSum %f'%ga4gh_attribute_value.double_value)
      if(attribute_name == 'PID'):
	assert len(ga4gh_attribute_values_list) == 1
	ga4gh_attribute_value = ga4gh_attribute_values_list[0]
	assert ga4gh_attribute_value.HasField('string_value')
	print('PID '+ga4gh_attribute_value.string_value)
      if(attribute_name == 'PL'):
	int_list = []
	for ga4gh_attribute_value in ga4gh_attribute_values_list:
	  assert ga4gh_attribute_value.HasField('int32_value')
	  int_list.append(ga4gh_attribute_value.int32_value)
	print('PL [ %s ]'%(','.join(str(x) for x in int_list)))
  print('END_NEW_CALL')


#for v in stream.parse(sys.stdin, variants_pb2.Variant):
  #print('Ha')

