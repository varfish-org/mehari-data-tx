#!/usr/bin/env bash

set -x

# Make the script directory available to all called scribes.
export SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Obtain genome release from environment variables.
export GENOME_RELEASE=${GENOME_RELEASE-grch37}

# Configure paths.
export DATA_DIR=${DATA_DIR-$HOME/mehari-data-tx}
export SEQREPO_ROOT_DIR=$DATA_DIR/seqrepo

# Parse configuration into environment variables.
export CDOT_RELEASE=$(jq -r ".$GENOME_RELEASE.cdot_release" config.json)
export CDOT_FILENAME=$(jq -r ".$GENOME_RELEASE.cdot_filename" config.json)
export ENSEMBL_RELEASE=$(jq -r ".$GENOME_RELEASE.ensembl_release" config.json)
export ENSEMBL_TOKEN=$(jq -r ".$GENOME_RELEASE.ensembl_token" config.json)

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

# Run the individual steps.
#bash -x $SCRIPT_DIR/download.sh
#bash -x $SCRIPT_DIR/seqrepo.sh
#bash -x $SCRIPT_DIR/pass-1.sh
bash -x $SCRIPT_DIR/fix-ncbi.sh
bash -x $SCRIPT_DIR/pass-2.sh
