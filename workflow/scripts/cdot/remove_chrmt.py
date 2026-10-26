#!/usr/bin/env python
"""Remove the chrMT transcripts of a RefSeq cdot JSON, as the Ensembl chrMT graft replaces them.

cdot names the RefSeq chrMT transcripts fake-rna-<gene>, without a sequence, and has no
chrMT tRNAs and rRNAs. The gene records with the NCBI Gene IDs of the grafted genes go too.
The build fails unless the graft has a transcript with the same gene, HGNC ID, exons and
strand for each removed transcript and the same HGNC ID for each removed gene record, and
unless no other transcript names a grafted gene. So RefSeq chrMT transcripts that the graft
lacks do not go unnoticed, and nothing off chrMT is removed. It also fails on a RefSeq chrMT
transcript with a real accession, as the graft may then no longer be needed.
"""

import gzip
import json
import sys
from contextlib import redirect_stderr

#: Contig name for chrMT
CHRMT_CONTIG = "NC_012920.1"


def location(tx: dict) -> tuple:
    """The contig, strand and exons of each genome build of a transcript."""
    return tuple(
        (build, alignment["contig"], alignment["strand"], tuple((exon[0], exon[1]) for exon in alignment["exons"]))
        for build, alignment in sorted(tx["genome_builds"].items())
    )


def describe(location: tuple) -> str:
    return ", ".join(
        f"{build} {contig}:{exons[0][0]}-{exons[-1][1]}{strand} ({len(exons)} exons)"
        for build, contig, strand, exons in location
    )


def main(cdot: dict, graft: dict) -> dict:
    grafted = {(tx.get("hgnc"), location(tx)): tx_id for tx_id, tx in graft["transcripts"].items()}
    removed = {
        tx_id: tx
        for tx_id, tx in cdot["transcripts"].items()
        if any(alignment["contig"] == CHRMT_CONTIG for alignment in tx["genome_builds"].values())
    }
    if accessions := sorted(tx_id for tx_id in removed if not tx_id.startswith("fake-rna-")):
        raise ValueError(f"RefSeq has chrMT transcripts with accessions {accessions}: revisit the Ensembl chrMT graft")
    unmatched = []
    for tx_id, tx in sorted(removed.items()):
        del cdot["transcripts"][tx_id]
        match = grafted.get((tx.get("hgnc"), location(tx)))
        print(
            f"Removed transcript {tx_id} ({tx.get('gene_name')}, gene {tx['gene_version']}, "
            f"HGNC:{tx.get('hgnc')}, {describe(location(tx))}), grafted: {match}",
            file=sys.stderr,
        )
        if match is None or graft["transcripts"][match]["gene_version"] != tx["gene_version"]:
            unmatched.append(tx_id)
    if unmatched:
        raise ValueError(f"The chrMT graft has no transcript with the gene, HGNC ID and location of {unmatched}")

    kept_genes = {tx["gene_version"] for tx in cdot["transcripts"].values()}
    if shared := kept_genes & graft["genes"].keys():
        raise ValueError(f"Transcripts off chrMT name the grafted chrMT genes {sorted(shared)}")
    removed_genes = sorted(cdot["genes"].keys() & graft["genes"].keys())
    for gene_id in removed_genes:
        gene = cdot["genes"].pop(gene_id)
        grafted_hgnc = graft["genes"][gene_id].get("hgnc")
        print(f"Removed gene {gene_id} ({gene.get('gene_symbol')}, HGNC:{gene.get('hgnc')})", file=sys.stderr)
        if gene.get("hgnc") != grafted_hgnc:
            raise ValueError(f"Gene {gene_id} has HGNC:{gene.get('hgnc')}, but the grafted gene has HGNC:{grafted_hgnc}")
    print(f"Removed {len(removed)} transcripts and {len(removed_genes)} genes", file=sys.stderr)
    return cdot


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    with gzip.open(snakemake.input.cdot, "rt") as file:
        cdot = json.load(file)
    with open(snakemake.input.cdot_mt, "rt") as file:
        graft = json.load(file)
    cdot = main(cdot, graft)
    with gzip.open(snakemake.output.cdot, "wt") as out:
        json.dump(cdot, out, indent=2)
