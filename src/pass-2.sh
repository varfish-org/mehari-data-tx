#!/usr/bin/env bash

# Also see comment at top of pass-1.sh.

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

set -x

mkdir -p $DATA_DIR/pass-2

sudo docker run \
    --volume "$DATA_DIR:$DATA_DIR:rw" \
    --volume "$(readlink -f $SCRIPT_DIR/..):$(readlink -f $SCRIPT_DIR/..):ro" \
    ghcr.io/varfish-org/mehari:${MEHARI_VERSION} \
    exec /usr/local/bin/mehari \
        db create \
        $(if [[ "$GENOME_RELEASE" == grch37 ]]; then \
            echo --path-mane-txs-tsv; \
            echo $DATA_DIR/tmp/mane-txs.tsv; \
        fi) \
        --path-out $DATA_DIR/pass-2/txs.bin.zst \
        --path-seqrepo-instance $DATA_DIR/seqrepo/master \
        --path-cdot-json $DATA_DIR/tmp/$GENOME_RELEASE/$CDOT_FILENAME \
        --path-cdot-json $SCRIPT_DIR/../data/cdot-0.2.24.ensembl.chrMT.$GENOME_RELEASE.gff3.json \
        --genome-release $GENOME_RELEASE
cd $DATA_DIR/pass-2

# Ensure that the output can be decompressed.
zstd -c -d txs.bin.zst > /dev/null

sha256sum txs.bin.zst > txs.bin.zst.sha256
sha256sum txs.bin.zst.report > txs.bin.zst.report.sha256
