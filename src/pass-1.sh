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

mehari db create txs \
    --path-out $DATA_DIR/pass-1/txs.bin.zst \
    --path-seqrepo-instance $DATA_DIR/seqrepo/master \
    --path-cdot-json $DATA_DIR/tmp/$GENOME_RELEASE/$CDOT_FILENAME \
    --genome-release $GENOME_RELEASE

# Ensure that the output can be decompressed.
zstd -c -d $DATA_DIR/pass-1/txs.bin.zst > /dev/null
