#!/bin/bash

# Define the description for the analysis
description="coloc credible set analysis April 1st with 25.03 platform version data"

# Define abbreviation for file name (simplified, no special characters)
abbreviation="2503_colocCredibleSet_DoE"

# Define the date for the log file
current_date=$(date '+%Y-%m-%d')

# Define the Python script to run
python_script="20250401_colocCredibleSet_DoE.py"

# Define the log file name (simplified)
log_filename="${current_date}_analysis_${abbreviation}.log"

# Write the description and other details to the log file
echo "--------------------------" > "$log_filename"
echo "Description: $description" >> "$log_filename"
echo "Abbreviation: $abbreviation" >> "$log_filename"
echo "Date: $current_date" >> "$log_filename"
echo "log file name: $log_filename" >> "$log_filename"
echo "--------------------------" >> "$log_filename"

# Run the runAnalysisGeneral.sh script and pass the description, log file name, and Python script as arguments
bash runAnalysisGeneral.sh "$description" "$current_date" "$python_script" "$log_filename" >> "$log_filename" 2>&1 &

# Tail the log file to see the last 1000 lines
tail -1000 "$log_filename"