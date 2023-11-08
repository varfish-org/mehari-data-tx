#!/usr/bin/env bash

set -x

# Download the GFF transcript models and FASTA sequences.

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

mkdir -p $DOWNLOAD_DIR
pushd $DOWNLOAD_DIR

urls=$(eval echo $(jq -r ".$GENOME_RELEASE | to_entries | map(select(.key | match(\"url\"))) | map(.value) | .[]" $SCRIPT_DIR/../config.json))

eval wget --continue --quiet $urls

popd
