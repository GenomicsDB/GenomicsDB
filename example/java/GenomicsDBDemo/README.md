Test GenomicsDB Demo Workspace

1. Set Environment Variable GENOMICSDB_DEMO_WS=/path/to/genomicsdb/demo/ws
2. Run `java -cp /path/to/genomicsdb-build/target/genomicsdb-1.5.0-SNAPSHOT-allinone-spark.jar /path/to/example/GenomicsTestGenomicsDBDemo.java`

e.g. Output
```
Summary Elapsed time: 568s
Interval: 2507743582-2507762731 calls=2962039

Summary for REF=="A": Elapsed time: 139s
Interval: 2507743582-2507762731 calls=400432

Summary for REF=="A" && ALT|="T": Elapsed time: 75s
Interval: 2507743582-2507762731 calls=82245

Summary for REF=="A" && ALT|="T" && GT&="1/1": Elapsed time: 98s
Interval: 2507743582-2507762731 calls=82245
```