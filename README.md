# mehari-data

Reproducible data builds for [mehari](https://github.com/bihealth/mehari).

This repository contains a snakemake workflow to build mehari transcript databases based on [SACGF/cdot](https://github.com/SACGF/cdot).

## Resulting Files

Databases are built for each combination of genome release (GRCh37, GRCh38) and reference source (ensembl, refseq), once for each source release listed in `config/config.yaml`:

- mehari-data `v0.14.0`
  - uses:
    - mehari `v0.47.0`
    - cdot: `v0.2.36`
  - `GRCh37-refseq-105.20220307`
      - genome release: GRCh37.p13
      - RefSeq annotation release: `105.20220307`
      - RefSeq assembly: [GCF\_000001405.25](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/)
      - resulting transcript database: `mehari-data-txs-GRCh37-refseq-105.20220307-0.14.0.bin.zst`
  - `GRCh37-ensembl-87`
      - genome release: GRCh37.p13
      - VEP/ENSEMBL release: `87`
      - RefSeq equivalent: [GCF\_000001405.25](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/)
      - resulting transcript database: `mehari-data-txs-GRCh37-ensembl-87-0.14.0.bin.zst`
  - `GRCh38-refseq-RS_2025_08`
      - genome release: GRCh38.p14
      - RefSeq annotation release: `RS_2025_08`
      - RefSeq assembly: [GCF\_000001405.40](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40/)
      - resulting transcript database: `mehari-data-txs-GRCh38-refseq-RS_2025_08-0.14.0.bin.zst`
  - `GRCh38-ensembl-115` and `GRCh38-ensembl-116`
      - genome release: GRCh38.p14
      - VEP/ENSEMBL release: `115` and `116`
      - RefSeq equivalent: [GCF\_000001405.40](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40/)
      - resulting transcript databases: `mehari-data-txs-GRCh38-ensembl-115-0.14.0.bin.zst` and `mehari-data-txs-GRCh38-ensembl-116-0.14.0.bin.zst`
  - `GRCh37-ensembl-87-and-refseq-105.20220307` and `GRCh38-ensembl-116-and-refseq-RS_2025_08`
      - the newest listed Ensembl and RefSeq databases of a genome release, merged
      - resulting transcript databases: `mehari-data-txs-GRCh37-ensembl-87-and-refseq-105.20220307-0.14.0.bin.zst` and `mehari-data-txs-GRCh38-ensembl-116-and-refseq-RS_2025_08-0.14.0.bin.zst`

Up to `v0.13.0`, the databases were named `mehari-data-txs-{genome release}-{ensembl,refseq,ensembl-and-refseq}-{version}.bin.zst`, with one release per source.


New builds of mehari-data will be considered when

- a new mehari version is released, OR
- a new cdot version is released, OR
- bugs are found and make a new release necessary.

Generally, only the data for the latest mehari protobuf schema is created.
