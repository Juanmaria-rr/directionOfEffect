#!/bin/bash

# Read the arguments passed from the calling script
description="$1"
current_date="$2"
python_script="$3"
log_filename="$4"

# Write the description and other details to the log file
echo "--------------------------" > "$log_filename"
echo "Description: $description" >> "$log_filename"
echo "Date: $current_date" >> "$log_filename"
echo "Python Script: $python_script" >> "$log_filename"
echo "--------------------------" >> "$log_filename"

# Run your Python script and log the output
python3 "$python_script" >> "$log_filename" 2>&1 &