#!/bin/bash

mkdir -p /tmp/vcfs
cp /opt/tests/inputs/vcfs/t0.vcf.gz /tmp/vcfs/
cp /opt/tests/inputs/vcfs/t0.vcf.gz.tbi /tmp/vcfs/
cp /opt/tests/inputs/vcfs/t1.vcf.gz /tmp/vcfs/
cp /opt/tests/inputs/vcfs/t1.vcf.gz.tbi /tmp/vcfs/
cp /opt/tests/inputs/vcfs/t2.vcf.gz /tmp/vcfs/
cp /opt/tests/inputs/vcfs/t2.vcf.gz.tbi /tmp/vcfs/
vcf2genomicsdb_init -w /tmp/ws -S /tmp/vcfs -n 0

if [ ! -f /tmp/ws/loader.json ]; then
  echo "vcf2genomicsdb_init failed"
  exit 1
fi

vcf2genomicsdb /tmp/ws/loader.json


tdb_file_count=`find /tmp/ws/ -name "*tdb"|wc -l`

if (( tdb_file_count < 10 )); then
  echo "vcf2genomicsdb failed, only found $tdb_file_count tdb files"
  exit 1
fi

tee -a /tmp/ws/query.json > /dev/null <<EOT
{
  "workspace" : "/tmp/ws",
  "array" : "allcontigs\$1\$3101976562",
  "query_column_ranges": [[[0,210306]]],
  "query_row_ranges": [[[0, 2]]],
  "query_attributes": ["GT","ALT","REF"],
  "reference_genome": "/opt/tests/inputs/chr1_10MB.fasta.gz",
  "callset_mapping_file": "/tmp/ws/callset.json",
  "vid_mapping_file": "/tmp/ws/vidmap.json",
  "produce_GT_field": true
}
EOT

variant_count=`gt_mpi_gather -j /tmp/ws/query.json --produce-Broad-GVCF| grep -v '#'| wc -l`

if (( variant_count != 4 )); then
  echo "gt_mpi_gather failed, expected 4 variants, found $variant_count"
  exit 1
fi
