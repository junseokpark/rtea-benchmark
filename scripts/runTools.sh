#!/bin/bash

# Simple tool execution script
# Configuration is loaded from the shared config.sh file

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the shared configuration file
if [ -f "${SCRIPT_DIR}/config.sh" ]; then
    source "${SCRIPT_DIR}/config.sh"
else
    echo "ERROR: Configuration file not found: ${SCRIPT_DIR}/config.sh"
    echo "Please create config.sh from config_template.sh and update the paths."
    exit 1
fi

module load singularity

# Execute JET

# Execute TEProf2

