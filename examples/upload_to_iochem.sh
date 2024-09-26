#!/bin/bash

# This script uploads Gaussian calculations to ioChem-BD. It has to be executed using `source upload_to_iochem.sh`.

# Configuration variables for this dataset
DIR_WITH_CALCS="/Volumes/Timo-3TB-WD-HDD/PhD/projects/DART/examples/Pd_Ni_Cross_Coupling/dev/output/gaussian_relaxed_complexes/batches"
IOCHEM_PROJECT="DART_Pd_Ni_intermediate"
# iochem shell configuration
IOCHEM_SHELL_PATH="$HOME/opt/ioChem-BD/shell/start-rep-shell"


# Helper function to extract the complex name from filename
get_complex_name_from_filename() {
    local filename=$(basename "$1")
    local name="${filename%_gaussian.log}"
    echo "$name"
}

# Function to upload gaussian calculations
upload_gaussian_calcs() {
    local log_file="$1"
    local com_file="${log_file%.log}.com"
    local complex_name=$(get_complex_name_from_filename "$log_file")
    local desc="$complex_name"
    local dir=$(dirname "$log_file")

    # Change to the directory with the log file
    cd "$dir"

    # Get only filename, otherwise iochem throws an error
    log_file=$(basename "$log_file")
    com_file=$(basename "$com_file")

#   echo "loadgauss -i \"$com_file\" -o \"$log_file\" -n \"$complex_name\" -d \"$desc\""  # For testing
    loadgauss -i "$com_file" -o "$log_file" -n "$complex_name" -d "$desc"
}

# Main execution starts here

# Load ioChem-BD API and navigate to the project
source $IOCHEM_SHELL_PATH
cdpro $IOCHEM_PROJECT

cwd=$(pwd)  # Save original directory, because the pwd will be changed

# Find all .log files and upload them
find "$DIR_WITH_CALCS" -type f -name '*.log' | while read -r log_file; do
    filename=$(basename "$log_file")
    echo "Uploading $filename"
    upload_gaussian_calcs "$log_file"
done

cd "$cwd" # Return to original directory