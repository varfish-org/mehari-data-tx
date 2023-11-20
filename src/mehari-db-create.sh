set -x

# Build Mehari database

# Check proper usage.
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 SOURCE RELEASE"
    exit 1
fi

SOURCE=$1
RELEASE=$2

# Make the script directory available to all called scripts.
export SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Get cdot version
export CDOT_VERSION=$(python $CDOT_DIR/generate_transcript_data/cdot_json.py --version)
# Get Mehari version
export MEHARI_VERSION=$(mehari --version | cut -d ' ' -f 2)
# Get label for version
export VERSION_LABEL=$(jq -r .version_label $SCRIPT_DIR/../config.json)

# Configure paths.
export DATA_DIR=${DATA_DIR-$HOME/mehari-data-tx}

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

# Create Mehari database.
mehari db create \
    --genome-release $RELEASE \
    --path-out $DATA_DIR/mehari-data-tx-$SOURCE-$RELEASE-$VERSION_LABEL+$MEHARI_VERSION+$CDOT_VERSION.bin.zst \
    $(if [[ "$SOURCE" == "refseq" ]] || [[ "$SOURCE" == "all" ]]; then \
        echo --path-cdot-json=$DATA_DIR/cdot-refseq-$VERSION_LABEL+$CDOT_VERSION.json.gz; \
    fi) \
    $(if [[ "$SOURCE" == "ensembl" ]] || [[ "$SOURCE" == "all" ]]; then \
        echo --path-cdot-json=$DATA_DIR/cdot-ensembl-$VERSION_LABEL+$CDOT_VERSION.json.gz; \
    fi) \
    $(if [[ "$RELEASE" == "grch37" ]]; then echo --path-mane-txs-tsv=$DATA_DIR/mane-txs.tsv; fi) \
    --path-seqrepo-instance $DATA_DIR/seqrepo/main

pushd $DATA_DIR
sha256sum mehari-data-tx-$SOURCE-$RELEASE-$VERSION_LABEL+$MEHARI_VERSION+$CDOT_VERSION.bin.zst \
> mehari-data-tx-$SOURCE-$RELEASE-$VERSION_LABEL+$MEHARI_VERSION+$CDOT_VERSION.bin.zst.sha256
popd
