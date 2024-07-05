#!/usr/bin/env python
import enum
import gzip
import json
import sys
from collections import defaultdict
from contextlib import redirect_stderr
from copy import deepcopy

import pandas as pd


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


LOCUS_TYPE_TO_BIOTYPE = {
    "complex locus constituent": "other",
    "endogenous retrovirus": "other",
    "fragile site": "other",
    "gene with protein product": "other",
    "immunoglobulin gene": "other",
    "immunoglobulin pseudogene": "other",
    "pseudogene": "pseudogene",
    "readthrough": "other",
    "region": "other",
    "RNA, cluster": "other",
    "RNA, long non-coding": "other",
    "RNA, micro": "other",
    "RNA, misc": "other",
    "RNA, ribosomal": "other",
    "RNA, small nuclear": "other",
    "RNA, small nucleolar": "other",
    "RNA, transfer": "other",
    "RNA, vault": "other",
    "RNA, Y": "other",
    "T cell receptor gene": "other",
    "T cell receptor pseudogene": "other",
    "unknown": "other",
    "virus integration site": "other",
}


def cdot_from_hgnc(cdot: str, hgnc: str):
    cdot = load_json(cdot)
    hgnc = load_json(hgnc)["response"]
    name_to_hgnc = defaultdict(lambda: defaultdict(str))
    id_to_hgnc = defaultdict(dict)
    hgnc_to_biotype = defaultdict(set)
    for record in hgnc["docs"]:
        hgnc_id = record["hgnc_id"][5:]  # strip "HGNC:" prefix
        for source in ["ensembl_gene_id", "symbol", "entrez_id"]:
            if r := record.get(source):
                name_to_hgnc[source][r] = hgnc_id
                if source == "entrez_id":
                    id_to_hgnc[r] = record
        if r := record.get("locus_type"):
            if b := LOCUS_TYPE_TO_BIOTYPE.get(r):
                if b != "other":
                    hgnc_to_biotype[hgnc_id].add(b)

    cdot_fixed = deepcopy(cdot)
    cdot_fixed = update_cdot(cdot, cdot_fixed, name_to_hgnc, hgnc_to_biotype)
    return cdot_fixed


def update_cdot(
    cdot: dict, update_target: dict, name_to_hgnc: dict, hgnc_to_biotype: dict
):
    """
    check transcript entries in cdot file for missing hgnc ids
    and add them from the hgnc file
    """
    report: list[tuple[str, str, str | None]] = []

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
            report.append(("gene", key, hgnc_id))
        biotype = gene.get("biotype", [])
        biotype = set(biotype) | hgnc_to_biotype.get(hgnc_id, set())
        gene["biotype"] = list(biotype)
        update_target["genes"].update({key: gene})

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
            report.append(("transcript", key, hgnc_id))
        biotype = tx.get("biotype", [])
        biotype = set(biotype) | hgnc_to_biotype.get(hgnc_id, set())
        tx["biotype"] = list(biotype)
        update_target["transcripts"].update({key: tx})

    return update_target, report


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    cdot, report = cdot_from_hgnc(
        snakemake.input.cdot,
        snakemake.input.hgnc,
    )

    with gzip.open(snakemake.output.cdot, "wt") as out:
        json.dump(cdot, out, indent=2)
        out.flush()
    pd.DataFrame(report, columns=["type", "key", "hgnc_id"]).to_csv(
        snakemake.output.report, sep="\t", index=False
    )
