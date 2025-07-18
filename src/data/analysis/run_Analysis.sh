#!/bin/bash

# Define the description for the analysis
description="L2G credible set analysis April 2nd with 25.03 platform version data"

# Define abbreviation for file name (simplified, no special characters)
abbreviation="2503_coherencyL2G"

# Define the date for the log file
current_date=$(date '+%Y-%m-%d')

# Define the Python script to run
python_script="20250402_newColoc_coherencyL2G.py"

# Define the log file name (simplified)
log_filename="${current_date}_analysis_${abbreviation}.log"

# Write the description and other details to the log file
echo "--------------------------" > "$log_filename"
echo "Description: $description" >> "$log_filename"
echo "Abbreviation: $abbreviation" >> "$log_filename"
echo "Date: $current_date" >> "$log_filename"
echo "Python Script: $python_script" >> "$log_filename"
echo "Log File: $log_filename" >> "$log_filename"
echo "--------------------------" >> "$log_filename"

# Run the Python script in the background and log output
nohup python3 "$python_script" >> "$log_filename" 2>&1 &

# Tail the log file in real-time
tail -1000 "$log_filename"