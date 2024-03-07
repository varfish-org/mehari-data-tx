# mehari-data

Reproducible data builds for [mehari](https://github.com/bihealth/mehari).

This repository contains scripts to build the transcript protobuf files for meharib based on [SACGF/cdot](https://github.com/SACGF/cdot).

## Resulting Files

The following explains the content and compatibility.

- mehari-data `v0.2.1`
  - compatible to: mehari `v0.4.1..`
  - `grch37` data file: `mehari-data-txs-grch37-0.2.1.bin.zst`
    - cdot: `v0.2.14`
    - genome release: GRCh37.p13
    - VEP/ENSEMBL equivalent: `r105`
    - RefSeq assembly: [GCF\_000001405.25](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/)
  - `grch38` data: `mehari-data-txs-grch38-0.2.1.bin.zst`
    - cdot: `v0.2.14`
    - genome release: GRCh38.p13
    - VEP/ENSEMBL equivalent: `r109`
    - RefSeq assembly: [GCF\_000001405.39](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/)

So, for example:

- `mehari-data-txs-grch37-0.2.1.bin.zst` is compatible to the mehari software `v0.4.1` and above;
- it was created from the cdot transcripts `v0.2.14`;
- these were built for GRCh37.p13 based on the VEP/ENSEMBL Release r105
- and the transcripts from the RefSeq assembly GCF\_000001405.25.

New builds of mehari-data will be considered when

- the new mehari software version changes protobuf schema, OR
- a new cdot has a new release corresponding to a VEP release, OR
- bugs are found and make a new release necessary.

Generally, only the data for the latest mehari protobuf schema is created.

## Utility Files

As RefSeq does not contain transcripts for mitochondrial genes, we graft the ENSEMBL transcripts over.

```
# python src/cdot_extract_chrmt.py \
      /tmp/cdot-0.2.24.ensembl.Homo_sapiens.GRCh37.87.gff3.json.gz \
  > data/cdot-0.2.24.ensembl.chrMT.grch37.gff3.json

# python src/cdot_extract_chrmt.py \
      /tmp/cdot-0.2.24.ensembl.Homo_sapiens.GRCh38.111.gff3.json.gz \
  > data/cdot-0.2.24.ensembl.chrMT.grch38.gff3.json
```
