#!/usr/bin/env bash

set -x

# Prepare the build by downloading the GFF and FASTA files, building the
# seqrepo, and converting the GFF files to cdot JSON files.

# Make the script directory available to all called scripts.
export SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Check proper usage.
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 GENOME_RELEASE"
    exit 1
fi

# Extract genome release argument to environment variable.
export GENOME_RELEASE=$1

# Get cdot version
export CDOT_VERSION=$(python $CDOT_DIR/generate_transcript_data/cdot_json.py --version)
# Get label for version
export VERSION_LABEL=$(jq -r .version_label $SCRIPT_DIR/../config.json)

# Configure paths.
export DATA_DIR=${DATA_DIR-$HOME/mehari-data-tx}
export DOWNLOAD_DIR=$DATA_DIR/downloads/$GENOME_RELEASE
export SEQREPO_ROOT_DIR=$DATA_DIR/seqrepo

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

# Cleanup $DATA_DIR if asked to do so.
if [ "${CLEANUP_DATA_DIR-false}" == "true" ]; then
    rm -rf $DATA_DIR
fi

# Run the individual steps.
bash $SCRIPT_DIR/download.sh
bash $SCRIPT_DIR/seqrepo.sh
# bash $SCRIPT_DIR/convert.sh
