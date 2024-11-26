# mehari-data

Reproducible data builds for [mehari](https://github.com/bihealth/mehari).

This repository contains a snakemake workflow to build mehari transcript databases based on [SACGF/cdot](https://github.com/SACGF/cdot).

## Resulting Files

Databases are built for each combination of genome release (GRCh37, GRCh38) and reference source (ensembl, refseq):

- mehari-data `v0.9.0`
  - uses:
    - mehari `v0.30.1`
    - cdot: `v0.2.27`
  - `GRCh37-refseq`
      - genome release: GRCh37.p13
      - VEP/ENSEMBL equivalent: `105`
      - RefSeq assembly: [GCF\_000001405.25](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/)
      - resulting transcript database: `mehari-0.26.1.GRCh37-refseq.txs.bin.zst`
  - `GRCh37-ensembl`
      - genome release: GRCh37.p13
      - VEP/ENSEMBL release: `105`
      - RefSeq equivalent: [GCF\_000001405.25](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/)
      - resulting transcript database: `mehari-0.26.1.GRCh37-ensembl.txs.bin.zst`
  - `GRCh38-refseq`
      - genome release: GRCh38.p13
      - VEP/ENSEMBL equivalent: `112`
      - RefSeq assembly: [GCF\_000001405.39](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/)
      - resulting transcript database: `mehari-0.26.1.GRCh38-refseq.txs.bin.zst`
  - `GRCh38-ensembl`
      - genome release: GRCh38.p13
      - VEP/ENSEMBL release: `112`
      - RefSeq equivalent: [GCF\_000001405.39](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/)
      - resulting transcript database: `mehari-0.26.1.GRCh38-ensembl.txs.bin.zst`


New builds of mehari-data will be considered when

- a new mehari version is released, OR
- a new cdot version is released, OR
- bugs are found and make a new release necessary.

Generally, only the data for the latest mehari protobuf schema is created.
