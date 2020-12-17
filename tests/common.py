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

import errno
import glob
import os
from os import environ
import subprocess
import sys
import zipfile
if sys.version_info[0] >= 3:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve

def __error(message):
    sys.stderr.write('Error: '+message+'\n')

def __error_exit(message):
    __error(message)
    sys.exit(-1)

def __find_genomicsdb_jar(target_dir, jar_file_name):
    jar_files=glob.glob(os.path.join(target_dir,jar_file_name))
    jar_file=''
    if (len(jar_files) == 1):
        jar_file=jar_files[0]
    else:
        if (len(jar_files) == 0):
            __error_exit('libraries matching '+jar_file_name+' not found')
        else:
            __error_exit('multiple libraries matching '+jar_file_name+ ' found')
    if (not os.path.isfile(jar_file)):
        __error_exit(jar_file+' found is not a library file')
    return jar_file

def setup_classpath(build_dir):
    target_dir=os.path.join(build_dir,'target')
    allinone_jar=__find_genomicsdb_jar(target_dir,'genomicsdb-*allinone.jar')
    examples_jar=__find_genomicsdb_jar(target_dir,'genomicsdb-*examples.jar')
    if 'CLASSPATH' in os.environ:
        classpath=os.environ['CLASSPATH']
    else:
        classpath=''
    if (len(classpath) > 0):
        environ["CLASSPATH"] = allinone_jar+os.pathsep+examples_jar
    else:
        environ["CLASSPATH"] = allinone_jar+os.pathsep+examples_jar+os.pathsep+classpath

def setup_jacoco(build_dir, build_type):
    jacoco= ''
    jacoco_report_cmd=''
    if (build_type.lower() == "coverage"):
        target_dir = os.path.join(build_dir, 'target')
        jacoco_cli_jar = os.path.join(build_dir,'lib/jacococli.jar')
        jacoco_agent_jar = os.path.join(build_dir,'lib/jacocoagent.jar')
        if (not os.path.isfile(jacoco_cli_jar) or not os.path.isfile(jacoco_agent_jar)):
            jacoco_zip_file = os.path.join(build_dir, "jacoco-0.8.5.zip")
            urlretrieve ("https://github.com/jacoco/jacoco/releases/download/v0.8.2/jacoco-0.8.2.zip", 
                                jacoco_zip_file)
            with zipfile.ZipFile(jacoco_zip_file, 'r') as zip_ref:
                zip_ref.extract('lib/jacocoagent.jar', build_dir)
                zip_ref.extract('lib/jacococli.jar', build_dir)
            os.remove(jacoco_zip_file)
        if (not os.path.isfile(jacoco_cli_jar) or not os.path.isfile(jacoco_agent_jar)):
            __error_exit('could not find jacoco cli or agent jar files')
        jacoco_reports_dir = os.path.join(target_dir, 'jacoco-reports')
        jacoco_dest_file = os.path.join(jacoco_reports_dir,"jacoco-ci.exec")
        jacoco = ' -javaagent:'+jacoco_agent_jar \
                 +'=destfile='+jacoco_dest_file \
                 +',includes=org.genomicsdb.*'
        try:
            os.makedirs(os.path.join(jacoco_reports_dir, 'jacoco-ci'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                __error_exit('could not create jacoco-reports dir:'+e.errno+' '+e.filename+' '+e.strerror)
        genomicsdb_classes_dir = os.path.join(target_dir, 'jacoco-classes')
        allinone_archive = zipfile.ZipFile(__find_genomicsdb_jar(target_dir,'genomicsdb-*allinone.jar'))
        for file in allinone_archive.namelist():
            if file.startswith('org/genomicsdb'):
                allinone_archive.extract(file, genomicsdb_classes_dir)
        jacoco_report_cmd = 'java -jar '+jacoco_cli_jar \
                            +' report '+jacoco_dest_file \
                            +' --classfiles '+genomicsdb_classes_dir \
                            +' --html '+jacoco_reports_dir \
                            +' --xml '+os.path.join(os.path.join(jacoco_reports_dir, 'jacoco-ci'),'jacoco-ci.xml')
    return jacoco, jacoco_report_cmd

def report_jacoco_coverage(jacoco_report_cmd):
    rc = -1
    if (len(jacoco_report_cmd) > 0):
        try:
            pid = subprocess.Popen(jacoco_report_cmd, shell=True, stdout=subprocess.PIPE);
            stdout_string = pid.communicate()[0]
            if(pid.returncode != 0):
                __error('jacoco report generation command:'+jacoco_report_cmd+' failed')
            else:
                rc = 0
        except Exception as e:
            __error('exception thrown while generating jacoco reports:'+str(e))
    else:
        rc = 0
    return rc
