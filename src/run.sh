#!/usr/bin/bash

export PYTHONPATH=${{ env.CDOT_DIR }}

# create gene information JSON file
EMAIL=manuel.holtgrewe@bih-charite.de \
bash ${{ env.CDOT_DIR }}/generate_transcript_data/gene_info.sh

# downlooad data, prepare cdot JSON files, and create seqrepos
bash src/prepare.sh grch38
bash src/prepare.sh grch37

# merge GRCh37 and GRCh38 files for refseq and ensembl
bash src/merge-cdot.sh ensembl
bash src/merge-cdot.sh refseq

# produce MANE transcript TSV file
bash src/mane-tx.sh

# build mehari data files
bash src/mehari-db-create.sh grch37
bash src/mehari-db-create.sh grch38
