#!/bin/sh

## do not forget to make "spark.stop" if a spark session have been previously created/started before running this analysis.
python3 run_geneEvidColocDoE_analysis.py 

#> run_test-$(date +"%d-%m-%Y").log 

# nohup bash 20240709_runAnalysis.sh > analysis_20240709.log 2>&1 &