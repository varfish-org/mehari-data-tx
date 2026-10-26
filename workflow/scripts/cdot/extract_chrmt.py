#!/usr/bin/env python
"""Helper script to alignments to chrMT from CDOT JSON.

This is used to graft ENSEMBL chrMT transcripts to the RefSeq transcripts
database. Their genes are keyed by NCBI Gene ID, like the other RefSeq genes.
"""

import gzip
import json
import sys
from contextlib import redirect_stderr


#: Contig name for chrMT
CHRMT_CONTIG = "NC_012920.1"


def key_genes_by_ncbi_gene_id(json_data: dict, hgnc_path: str):
    """Key the genes by the NCBI Gene ID that HGNC names for their HGNC ID."""
    with open(hgnc_path, "rt") as inputf:
        docs = json.load(inputf)["response"]["docs"]
    ncbi_gene_ids = {doc["hgnc_id"].removeprefix("HGNC:"): doc.get("entrez_id") for doc in docs}
    new_gene_ids = {}
    for gene_id, gene in json_data["genes"].items():
        ncbi_gene_id = ncbi_gene_ids.get(gene.get("hgnc"))
        if not ncbi_gene_id:
            raise ValueError(
                f"HGNC names no NCBI Gene ID for {gene_id} ({gene.get('gene_symbol')}, HGNC:{gene.get('hgnc')})"
            )
        if ncbi_gene_id in new_gene_ids.values():
            raise ValueError(f"HGNC names the NCBI Gene ID {ncbi_gene_id} of another grafted gene for {gene_id}")
        print(f"Keyed gene {gene_id} ({gene.get('gene_symbol')}) by NCBI Gene ID {ncbi_gene_id}", file=sys.stderr)
        new_gene_ids[gene_id] = ncbi_gene_id
    json_data["genes"] = {new_gene_ids[gene_id]: gene for gene_id, gene in json_data["genes"].items()}
    for tx in json_data["transcripts"].values():
        tx["gene_version"] = new_gene_ids[tx["gene_version"]]


def main(json_path: str, hgnc_path: str, file=sys.stdout):
    #: Transcript IDs (with version) to keep.
    tx_ids: set[str] = set()
    #: Gene IDs (with version) tok eep.
    gene_ids: set[str] = set()

    print(f"processing {json_path}...", file=sys.stderr)

    print(f"- loading {json_path}", file=sys.stderr)
    if json_path.endswith(".gz"):
        with gzip.open(json_path, "rt") as inputf:
            json_data = json.load(inputf)
    else:
        with open(json_path, "rt") as inputf:
            json_data = json.load(inputf)

    print("- processing...", file=sys.stderr)
    for idx, tx in enumerate(json_data["transcripts"].values()):
        if idx > 0 and idx % 1000 == 0:
            print(f"  ... processed {idx} transcripts", file=sys.stderr)
        for genome_build in tx["genome_builds"].values():
            if genome_build["contig"] == CHRMT_CONTIG:
                # print(json.dump(tx, sys.stderr, indent=2))
                # print("\n", file=sys.stderr)
                tx_ids.add(tx["id"])
                gene_ids.add(tx["gene_version"])
    print(f"... done processing {json_path}", file=sys.stderr)

    print("reducing to chrMT...", file=sys.stderr)
    for gene_id in list(json_data["genes"].keys()):
        if gene_id not in gene_ids:
            json_data["genes"].pop(gene_id)
    for tx_id in list(json_data["transcripts"].keys()):
        if tx_id not in tx_ids:
            json_data["transcripts"].pop(tx_id)
        else:
            tx = json_data["transcripts"][tx_id]
            for genome_build in tx["genome_builds"].values():
                if tags := genome_build.get("tag"):
                    tags = tags.split(",")
                else:
                    tags = []
                genome_build["tag"] = ",".join([*tags, "EnsemblGraft"])
                print("Updated transcript {tx_id} tags to {tags}".format(tx_id=tx_id, tags=genome_build["tag"]), file=sys.stderr)
    print("... done reducing to chrMT", file=sys.stderr)

    key_genes_by_ncbi_gene_id(json_data, hgnc_path)

    json.dump(json_data, file, indent=2)
    print(file=file, flush=True)


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    with open(snakemake.output.cdot, "w") as out:
        main(json_path=snakemake.input.cdot, hgnc_path=snakemake.input.hgnc, file=out)
