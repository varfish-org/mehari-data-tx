#!/usr/bin/env bash

# Build mehari transcript protobuf file (pass 1).
#
# Because the NCBI transcripts are not stable, fetching the transcript sequences
# will fail for certain transcripts.  This will be apparent from the report file
# written out by this pass.  The `fix-seqrepo.sh` script will attempt to work
# around this issue by fetching the sequences directly from NCBI and then `pass-2.sh`
# will work without such failures.

set -euo pipefail

export USE_CONDA_ENV=${USE_CONDA_ENV-mehari-data}
export DATA_DIR=${DATA_DIR-/tmp/mehari-data}
export SEQREPO_ROOT_DIR=$DATA_DIR/seqrepo
export PATH=$PATH:/home/runner/micromamba-root/envs/mehari-data/bin

if [[ ! -z "${ACTIVATE_CONDA-}" ]]; then
    source $ACTIVATE_CONDA $USE_CONDA_ENV
fi

set -x

mkdir -p $DATA_DIR/pass-1

if [[ "$GENOME_RELEASE" == grch37 ]]; then
    python src/cdot_json_to_tags.py \
        $DATA_DIR/tmp/$GENOME_RELEASE/$CDOT_FILENAME_38 \
    > $DATA_DIR/tmp/mane-txs.tsv
fi

mehari db create \
    $(if [[ "$GENOME_RELEASE" == grch37 ]]; then \
        echo --path-mane-txs-tsv; \
        echo $DATA_DIR/tmp/mane-txs.tsv; \
    fi) \
    --path-out $DATA_DIR/pass-1/txs.bin.zst \
    --path-seqrepo-instance $DATA_DIR/seqrepo/master \
    --path-cdot-json $DATA_DIR/tmp/$GENOME_RELEASE/$CDOT_FILENAME \
    --path-cdot-json $SCRIPT_DIR/../data/cdot-0.2.23.ensembl.chrMT.$GENOME_RELEASE.gff3.json \
    --genome-release $GENOME_RELEASE

# Ensure that the output can be decompressed.
zstd -c -d $DATA_DIR/pass-1/txs.bin.zst > /dev/null
