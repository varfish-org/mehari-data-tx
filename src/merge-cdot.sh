set -x

# Merge the cdot files for the given source.

# Make the script directory available to all called scripts.
export SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Check proper usage.
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 SOURCE"
    exit 1
fi

# Extract source argument to environment variable.
export DATA_SOURCE=$1

# Get cdot version
export CDOT_VERSION=$(python $CDOT_DIR/generate_transcript_data/cdot_json.py --version)
# Get label for version
export VERSION_LABEL=$(jq -r .version_label $SCRIPT_DIR/../config.json)

# Configure paths.
export DATA_DIR=${DATA_DIR-$HOME/mehari-data-tx}

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

python $CDOT_DIR/generate_transcript_data/cdot_json.py combine_builds \
    --grch37 $DATA_DIR/cdot-$DATA_SOURCE-grch37-$VERSION_LABEL-$CDOT_VERSION.json.gz \
    --grch38 $DATA_DIR/cdot-$DATA_SOURCE-grch38-$VERSION_LABEL-$CDOT_VERSION.json.gz \
    --output $DATA_DIR/cdot-$DATA_SOURCE-$VERSION_LABEL-$CDOT_VERSION.json.gz
