#!/usr/bin/env python
"""
This is used to graft ENSEMBL transcripts to the RefSeq transcripts
database.
"""

import gzip
import json
import sys

import pandas as pd
from contextlib import redirect_stderr
from snakemake.script import snakemake

REFSEQ_KEY = "RefSeq mRNA ID"
ENSEMBL_GENE_KEY = "Gene stable ID version"
ENSEMBL_TRANSCRIPT_KEY = "Transcript stable ID version"


def main(cdot, lookup):
    ensembl_genes = set(lookup[ENSEMBL_GENE_KEY])
    ensembl_transcripts = set(lookup[ENSEMBL_TRANSCRIPT_KEY])
    for gene_id in list(cdot["genes"].keys()):
        if gene_id not in ensembl_genes:
            cdot["genes"].pop(gene_id)
    for tx_id in list(cdot["transcripts"].keys()):
        if tx_id not in ensembl_transcripts:
            cdot["transcripts"].pop(tx_id)
    return cdot


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    with gzip.open(snakemake.input.cdot, "r") as f:
        cdot = json.load(f)
    lookup = pd.read_csv(snakemake.input.lookup, sep="\t")
    if not lookup.empty:
        lookup = lookup.set_index(REFSEQ_KEY, drop=False)
        cdot = main(cdot, lookup)
    else:
        print("Empty lookup table, skipping grafting.", file=sys.stderr)

    with gzip.open(snakemake.output.cdot, "wt") as out:
        json.dump(cdot, out, indent=2)
