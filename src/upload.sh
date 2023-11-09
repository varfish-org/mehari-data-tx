#!/usr/bin/bash

set -x

# Make the script directory available to all called scripts.
export SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Get cdot version
export CDOT_VERSION=$(python $CDOT_DIR/generate_transcript_data/cdot_json.py --version)
# Get label for version
export VERSION_LABEL=$(jq -r .version_label $SCRIPT_DIR/../config.json)

# Configure paths.
export DATA_DIR=${DATA_DIR-$HOME/mehari-data-tx}

# Upload the data...
gh release upload \
    $RELEASE_NAME \
    $DATA_DIR/cdot-ensembl-grch37-$VERSION_LABEL-$CDOT_VERSION.json.gz{,.sha256} \
    $DATA_DIR/cdot-refseq-grch37-$VERSION_LABEL-$CDOT_VERSION.json.gz{,.sha256} \
    $DATA_DIR/cdot-ensembl-grch38-$VERSION_LABEL-$CDOT_VERSION.json.gz{,.sha256} \
    $DATA_DIR/cdot-refseq-grch38-$VERSION_LABEL-$CDOT_VERSION.json.gz{,.sha256} \
    $DATA_DIR/mane-txs.tsv{,.sha256} \
    $DATA_DIR/mehari-data-tx-*.bin.zstd*
