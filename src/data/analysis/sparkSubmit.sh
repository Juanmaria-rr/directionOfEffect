#!/bin/bash

# Define a description for the analysis
#description="save matrices with coincidences and discrepancies, latest"
description="L2G benchmark with 25.03 platform version data // coloc from 25.03 release // included resolving discrepancies with N "

# Define an abbreviation for the file name (simplified, no special characters)
#abbreviation="matricesDiscr"
abbreviation="L2Gbenchmark_2503rel"

# Define the current date for the log file
current_date=$(date '+%Y-%m-%d')
#current_date=$(date '+%Y-%m-%d_%H-%M-%S')

# Define the Python script to run
#python_script="saveMatrices.py"
python_script="20250507_newColoc_coherencyL2G.py"

# Define the log file name
log_filename="${current_date}_analysis_${abbreviation}.log"

# Write the description and other details to the log file
echo "--------------------------" > "$log_filename"
echo "Description: $description" >> "$log_filename"
echo "Abbreviation: $abbreviation" >> "$log_filename"
echo "Date: $current_date" >> "$log_filename"
echo "Python Script: $python_script" >> "$log_filename"
echo "Log File Name: $log_filename" >> "$log_filename"
echo "--------------------------" >> "$log_filename"

# Run spark-submit and log the output
spark-submit "$python_script" >> "$log_filename" 2>&1 &

# Tail the log file to see the output in real-time
echo "tail -10 $log_filename"
#tail -f "$log_filename"
