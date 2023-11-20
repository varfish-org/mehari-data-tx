#!/usr/bin/bash

export PYTHONPATH=$CDOT_DIR

# create gene information JSON file
EMAIL=manuel.holtgrewe@bih-charite.de \
bash $CDOT_DIR/generate_transcript_data/gene_info.sh

# downlooad data, prepare cdot JSON files, and create seqrepos
bash src/prepare.sh grch38
bash src/prepare.sh grch37

# produce MANE transcript TSV file
bash src/mane-tx.sh

# merge GRCh37 and GRCh38 files for refseq and ensembl
bash src/merge-cdot.sh ensembl
bash src/merge-cdot.sh refseq

# build mehari data files
for source in ensembl refseq all; do
    for release in grch37 grch38; do
        bash src/mehari-db-create.sh $source
        bash src/mehari-db-create.sh $release
    done
done
