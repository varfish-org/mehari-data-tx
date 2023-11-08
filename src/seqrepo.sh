#!/usr/bin/env bash

set -x

# Create a seqrepo with the downloaded FASTA files.

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

if [[ ! -e $SEQREPO_ROOT_DIR/main ]]; then
    mkdir -p $SEQREPO_ROOT_DIR
    seqrepo init -i main
fi

# seqrepo load is too verbose and we cannot silence it
2>&1 seqrepo \
    load --instance-name main --namespace ENSEMBL \
    $(find $DOWNLOAD_DIR -name '*.fna.gz' -or -name '*.fa.gz' | grep -v _rna) \
| tail

# seqrepo load is too verbose and we cannot silence it
2>&1 seqrepo \
    load --instance-name main --namespace RefSeq \
    $DOWNLOAD_DIR/*_rna.fna.gz \
| tail
