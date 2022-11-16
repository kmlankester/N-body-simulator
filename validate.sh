#!/bin/bash
set -e

# Check if the needed files are in the repository
for FILE in step-1.cpp step-2.cpp step-3.cpp step-4.cpp report.pdf
do
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        echo "ERROR:" $FILE " does not exist."
        exit 1
    fi
done
