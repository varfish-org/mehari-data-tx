#!/usr/bin/env bash

set -x

# Convert GFF transcript models to cdot JSON files.

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

ensembl_out=$DATA_DIR/cdot-ensembl-$GENOME_RELEASE-$VERSION_LABEL-$CDOT_VERSION.json.gz
refseq_out=$DATA_DIR/cdot-refseq-$GENOME_RELEASE-$VERSION_LABEL-$CDOT_VERSION.json.gz

if [[ ! -e $ensembl_out ]]; then
    ensembl_url=$(jq -r ".$GENOME_RELEASE.ensembl_url_gff" $SCRIPT_DIR/../config.json)
    ensembl_file=$(echo $ensembl_url | rev | cut -d / -f 1 | rev)
    python $CDOT_DIR/generate_transcript_data/cdot_json.py gff3_to_json \
        --url=$ensembl_url \
        --genome-build=$(echo $GENOME_RELEASE | sed -e 's/grch37/GRCh37/' -e 's/grch38/GRCh38/') \
        --gene-info-json=$SCRIPT_DIR/../Homo_sapiens.gene-info-$CDOT_VERSION.json.gz \
        --output=$ensembl_out \
        $DOWNLOAD_DIR/$ensembl_file
    pushd $DATA_DIR
    sha256sum $(basename $ensembl_out) >$(basename $ensembl_out).sha256
    popd
fi

if [[ ! -e $refseq_out ]]; then
    refseq_url=$(jq -r ".$GENOME_RELEASE.refseq_url_gff" $SCRIPT_DIR/../config.json)
    refseq_file=$(echo $refseq_url | rev | cut -d / -f 1 | rev)
    python $CDOT_DIR/generate_transcript_data/cdot_json.py gff3_to_json \
        --url=$refseq_url \
        --genome-build=$(echo $GENOME_RELEASE | sed -e 's/grch37/GRCh37/' -e 's/grch38/GRCh38/') \
        --gene-info-json=$SCRIPT_DIR/../Homo_sapiens.gene-info-$CDOT_VERSION.json.gz \
        --output=$refseq_out \
        $DOWNLOAD_DIR/$refseq_file
    pushd $DATA_DIR
    sha256sum $(basename $refseq_out) >$(basename $refseq_out).sha256
    popd
fi
