set -x

# Create MANE transcript TSV files.

# Make the script directory available to all called scripts.
export SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Get cdot version
export CDOT_VERSION=$(python $CDOT_DIR/generate_transcript_data/cdot_json.py --version)
# Get label for version
export VERSION_LABEL=$(jq -r .version_label $SCRIPT_DIR/../config.json)

# Configure paths.
export DATA_DIR=${DATA_DIR-$HOME/mehari-data-tx}

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

python src/cdot_json_to_tags.py \
    $DATA_DIR/cdot-ensembl-grch38-$VERSION_LABEL-$CDOT_VERSION.json.gz \
    $DATA_DIR/cdot-refseq-grch38-$VERSION_LABEL-$CDOT_VERSION.json.gz \
> $DATA_DIR/mane-txs.tsv

pushd $DATA_DIR
sha256sum mane-txs.tsv >mane-txs.tsv.sha256
popd

