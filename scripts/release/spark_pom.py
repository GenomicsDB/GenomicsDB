# The MIT License (MIT)
# Copyright (c) 2026 Tim Fennell
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""Writes the pom of org.genomicsdb:genomicsdb-spark to stdout, made from genomicsdb's released pom given as the
argument: the same project description, with genomicsdb-spark's own coordinates and a single dependency, on genomicsdb
at the same version, which brings in the rest.

usage: spark_pom.py <genomicsdb-<version>.pom>
"""
import sys
import xml.etree.ElementTree as ET

NAMESPACE = "http://maven.apache.org/POM/4.0.0"
# The elements of genomicsdb's pom that describe the project, as Maven Central requires, rather than how it is built
PROJECT_ELEMENTS = {"modelVersion", "groupId", "artifactId", "version", "name", "description", "url", "licenses",
                    "developers", "scm"}


def element(tag):
    return f"{{{NAMESPACE}}}{tag}"


ET.register_namespace("", NAMESPACE)
ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")
project = ET.parse(sys.argv[1]).getroot()
for child in list(project):
    if child.tag.removeprefix(element("")) not in PROJECT_ELEMENTS:
        project.remove(child)

group_id = project.find(element("groupId")).text
version = project.find(element("version")).text
project.find(element("artifactId")).text = "genomicsdb-spark"
project.find(element("name")).text = f"{group_id}:genomicsdb-spark"
project.find(element("description")).text = "Spark API for GenomicsDB"
dependency = ET.SubElement(ET.SubElement(project, element("dependencies")), element("dependency"))
for tag, text in (("groupId", group_id), ("artifactId", "genomicsdb"), ("version", version)):
    ET.SubElement(dependency, element(tag)).text = text

ET.indent(project)
ET.ElementTree(project).write(sys.stdout.buffer, encoding="UTF-8", xml_declaration=True)
sys.stdout.buffer.write(b"\n")
