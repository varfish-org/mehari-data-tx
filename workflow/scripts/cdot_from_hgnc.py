#!/usr/bin/env python
"""Helper script to extract transcript identifier to tags from GRCh38.

The resulting file will be a TSV file and have the columns (1) tx_id,
(2) gene symbol, and (3) tags.  This can be used as the input for mehari
database building for applying the tags.
"""

import gzip
import json
import sys
from collections import defaultdict
from contextlib import redirect_stderr


def load_json(path: str) -> dict:
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as inputf:
            return json.load(inputf)
    else:
        with open(path, "rt") as inputf:
            return json.load(inputf)


def update_cdot_hgnc_ids(cdot: str, hgnc: str, file=sys.stdout):
    cdot = load_json(cdot)
    hgnc = load_json(hgnc)["response"]
    name_to_hgnc = defaultdict(lambda: defaultdict(str))
    id_to_hgnc = defaultdict(dict)
    for record in hgnc["docs"]:
        hgnc_id = record["hgnc_id"][5:]  # strip "HGNC:" prefix
        for source in ["ensembl_gene_id", "symbol", "entrez_id"]:
            if r := record.get(source):
                name_to_hgnc[source][r] = hgnc_id[5:]
                if source == "entrez_id":
                    id_to_hgnc[r] = record
    hgnc_refseq_ids = set(id_to_hgnc.keys())
    cdot_refseq_ids = set(cdot["genes"].keys())

    hgnc_cdot = {
        "cdot_version": cdot["cdot_version"],
        "genome_builds": cdot["genome_builds"],
        "genes": {},
        "transcripts": {},
    }

    # Generate cdot records from hgnc information for records missing from the given cdot file
    for refseq_id in hgnc_refseq_ids - cdot_refseq_ids:
        r = id_to_hgnc[refseq_id]
        record = {
            "aliases": ", ".join(
                s
                for s in (
                    r.get("symbol", ""),
                    ", ".join(r.get("alias_symbol", [])),
                    r.get("ensembl_gene_id", ""),
                    r.get("entrez_id", ""),
                    ", ".join(r.get("refseq_accession", [])),
                )
                if s
            ),
            "biotype": [],  # r.get("locus_type", "")
            "description": r.get("name", ""),
            "gene_symbol": r.get("symbol", ""),
            "hgnc": r.get("hgnc_id", "HGNC:")[5:],
            "map_location": r.get("location", ""),
            "summary": "missing from cdot; added from hgnc",
            "url": "https://example.com",
        }
        hgnc_cdot["genes"][refseq_id] = record

    # check transcript entries in cdot file for missing hgnc ids
    # and add them from the hgnc file
    for key, tx in cdot["transcripts"].items():
        for genome_build in tx.get("genome_builds", dict()).values():
            if genome_build.get("tag"):
                gene = tx["gene_name"]
                hgnc_id = tx.get("hgnc", None)
                if not hgnc_id:
                    print(f"Transcript without HGNC ID:", gene, file=sys.stderr)
                    for source in ["entrez_id", "ensembl_gene_id", "symbol"]:
                        if hgnc_id := name_to_hgnc[source].get(gene, None):
                            print("→ Corresponding HGNC ID:", hgnc_id, file=sys.stderr)
                            tx["hgnc"] = hgnc_id
                            hgnc_cdot["transcripts"].update({key: tx})
                            break

    json.dump(hgnc_cdot, file, indent=2)


with gzip.open(snakemake.output.cdot, "wt") as out:
    with open(snakemake.log[0], "w") as log, redirect_stderr(log):
        update_cdot_hgnc_ids(snakemake.input.cdot, snakemake.input.hgnc, file=out)
