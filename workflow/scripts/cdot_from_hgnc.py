#!/usr/bin/env python
import enum
import gzip
import json
import sys
from collections import defaultdict
from contextlib import redirect_stderr
from copy import deepcopy


def load_json(path: str) -> dict:
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as inputf:
            return json.load(inputf)
    else:
        with open(path, "rt") as inputf:
            return json.load(inputf)


class Mode(enum.StrEnum):
    create = "create"
    update = "update"


def cdot_from_hgnc(cdot: str, hgnc: str):
    cdot = load_json(cdot)
    hgnc = load_json(hgnc)["response"]
    name_to_hgnc = defaultdict(lambda: defaultdict(str))
    id_to_hgnc = defaultdict(dict)
    for record in hgnc["docs"]:
        hgnc_id = record["hgnc_id"][5:]  # strip "HGNC:" prefix
        for source in ["ensembl_gene_id", "symbol", "entrez_id"]:
            if r := record.get(source):
                name_to_hgnc[source][r] = hgnc_id
                if source == "entrez_id":
                    id_to_hgnc[r] = record

    cdot_fixed = deepcopy(cdot)
    cdot_fixed = update_cdot(cdot, cdot_fixed, name_to_hgnc)
    return cdot_fixed


def update_cdot(cdot: dict, update_target: dict, name_to_hgnc: dict):
    """
    check transcript entries in cdot file for missing hgnc ids
    and add them from the hgnc file
    """

    for key, gene in cdot["genes"].items():
        hgnc_id = gene.get("hgnc", None)
        if not hgnc_id:
            print("Gene without HGNC ID:", key, file=sys.stderr)
            for source in ["entrez_id", "ensembl_gene_id", "symbol"]:
                if hgnc_id := name_to_hgnc[source].get(
                    key, name_to_hgnc[source].get(gene["gene_symbol"], None)
                ):
                    print("→ Corresponding HGNC ID:", hgnc_id, file=sys.stderr)
                    gene["hgnc"] = hgnc_id
                    update_target["genes"].update({key: gene})
                    break

    for key, tx in cdot["transcripts"].items():
        gene = tx["gene_name"]
        hgnc_id = tx.get("hgnc", None)
        if not hgnc_id:
            print("Transcript without HGNC ID:", gene, file=sys.stderr)
            for source in ["entrez_id", "ensembl_gene_id", "symbol"]:
                if hgnc_id := name_to_hgnc[source].get(gene, None):
                    print("→ Corresponding HGNC ID:", hgnc_id, file=sys.stderr)
                    tx["hgnc"] = hgnc_id
                    update_target["transcripts"].update({key: tx})
                    break

    return update_target


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    cdot = cdot_from_hgnc(
        snakemake.input.cdot,
        snakemake.input.hgnc,
    )

    with gzip.open(snakemake.output.cdot, "wt") as out:
        json.dump(cdot, out, indent=2)
        out.flush()
