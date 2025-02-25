#!/usr/bin/env python
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
    id_to_hgnc_record = defaultdict(dict)
    hgnc_to_biotype = defaultdict(set)
    hgnc_to_symbol = defaultdict(str)
    sources = [
        "entrez_id",
        "ensembl_gene_id",
        "refseq_accession",
        "mane_select",
    ]
    # gene name to hgnc id
    # order is important, as the first match is taken
    symbols = ["symbol", "alias_symbol", "gene_name", "gene_symbol"]

    sources += symbols
    for record in hgnc["docs"]:
        hgnc_id = record["hgnc_id"][5:]  # strip "HGNC:" prefix
        id_to_hgnc_record[hgnc_id] = record
        for source in sources:
            if r := record.get(source):
                if source == "ensembl_gene_id" and "." in r:
                    acc = r.rsplit(".", 1)[0]
                    name_to_hgnc[source][acc] = hgnc_id
                if source == "mane_select":
                    ensembl = list(filter(lambda x: x.startswith("ENST"), r)) or []
                    refseq = list(filter(lambda x: x.startswith("NM_"), r)) or []
                    if ensembl:
                        for entry in ensembl:
                            name_to_hgnc["mane_select_ensembl"][entry] = hgnc_id
                    if refseq:
                        for entry in refseq:
                            name_to_hgnc["mane_select_refseq"][entry] = hgnc_id
                else:
                    if not isinstance(r, list):
                        r = [r]
                    for r_ in r:
                        name_to_hgnc[source][r_] = hgnc_id
                if source in {"symbol", "alias_symbol", "gene_symbol", "gene_name"}:
                    if not hgnc_to_symbol[hgnc_id]:
                        hgnc_to_symbol[hgnc_id] = r

        if r := record.get("locus_type"):
            if b := LOCUS_TYPE_TO_BIOTYPE.get(r):
                if b != "other":
                    hgnc_to_biotype[hgnc_id].add(b)

    cdot_fixed = deepcopy(cdot)
    cdot_fixed = update_cdot(
        cdot,
        cdot_fixed,
        name_to_hgnc,
        hgnc_to_biotype,
        hgnc_to_symbol,
        id_to_hgnc_record,
    )
    return cdot_fixed


def update_cdot(
    cdot: dict,
    update_target: dict,
    name_to_hgnc: dict,
    hgnc_to_biotype: dict,
    hgnc_to_symbol: dict,
    id_to_hgnc_record: dict,
):
    """
    check transcript entries in cdot file for missing hgnc ids
    and add them from the hgnc file
    """
    report: list[tuple[str, str, str, str | None, str]] = []

    for key, gene in cdot["genes"].items():
        _hgnc_id = None
        gene_hgnc_id = gene.get("hgnc", None)
        if "." in key:
            acc = key.rsplit(".", 1)[0]
        else:
            acc = key
        if not gene_hgnc_id:
            print("Gene without HGNC ID:", key, file=sys.stderr)
            for source in ["entrez_id", "ensembl_gene_id", "symbol"]:
                if hgnc_id := name_to_hgnc[source].get(
                    key, name_to_hgnc[source].get(acc, None)
                ):
                    print(f"→ Corresponding HGNC ID (via {source}):", hgnc_id, file=sys.stderr)
                    gene["hgnc"] = hgnc_id
                    _hgnc_id = hgnc_id
                    update_target["genes"].update({key: gene})
                    break
            report.append(("gene", "hgnc", key, None, _hgnc_id))
        # else:
        #     for source in ["entrez_id", "ensembl_gene_id", "symbol"]:
        #         if hgnc_id := name_to_hgnc[source].get(
        #             key, name_to_hgnc[source].get(acc, None)
        #         ):
        #             if hgnc_id == gene_hgnc_id:
        #                 continue
        #             print("Gene with outdated HGNC ID:", key, file=sys.stderr)
        #             print(
        #                 f"→ Updated HGNC ID from {gene_hgnc_id} to {hgnc_id}",
        #                 file=sys.stderr,
        #             )
        #             gene["hgnc"] = hgnc_id
        #             update_target["genes"].update({key: gene})
        #             report.append(("gene", "hgnc", key, gene_hgnc_id, hgnc_id))
        #             break
        gene_symbol = gene.get("gene_symbol", None)
        if not gene_symbol:
            print("Gene without gene symbol:", key, file=sys.stderr)
            if symbol := hgnc_to_symbol.get(_hgnc_id, None):
                symbol = symbol[0]
                print("→ Corresponding gene symbol:", symbol, file=sys.stderr)
                gene["gene_symbol"] = symbol
                update_target["genes"].update({key: gene})
            report.append(("gene", "gene_symbol", key, None, symbol))

        biotype_orig = list(sorted(gene.get("biotype", [])))
        biotype = list(sorted(set(biotype_orig) | hgnc_to_biotype.get(_hgnc_id, set())))
        if biotype_orig != biotype:
            gene["biotype"] = biotype
            update_target["genes"].update({key: gene})
            report.append(
                ("gene", "biotype", key, ",".join(biotype_orig), ",".join(biotype))
            )

    for key, tx in cdot["transcripts"].items():
        keys = (key, tx["id"], tx["gene_name"])
        if "." in key:
            keys += (key.rsplit(".", 1)[0],)
        transcript_hgnc_id = tx.get("hgnc", None)
        sources = [
            "entrez_id",
            "ensembl_gene_id",
            "refseq_accession",
            "mane_select_ensembl",
            "mane_select_refseq",
        ]
        symbols = ["symbol", "alias_symbol", "gene_name", "gene_symbol"]
        sources += symbols
        _hgnc_id = transcript_hgnc_id
        if not transcript_hgnc_id:
            print("Transcript without HGNC ID:", keys, file=sys.stderr)
            # order matters -- all matches are reported but only the last match is kept
            for source in reversed(sources):
                for k in keys:
                    if hgnc_id := name_to_hgnc[source].get(k, None):
                        print("→ Corresponding HGNC ID:", hgnc_id, file=sys.stderr)
                        tx["hgnc"] = hgnc_id
                        _hgnc_id = hgnc_id
                        update_target["transcripts"].update({key: tx})
            report.append(("transcript", "hgnc", key, None, _hgnc_id))
        # else:
        #     for source in sources:
        #         for k in keys:
        #             if hgnc_id := name_to_hgnc[source].get(k, None):
        #                 if hgnc_id == transcript_hgnc_id:
        #                     continue
        #                 print("Transcript with outdated HGNC ID:", k, file=sys.stderr)
        #                 print(
        #                     f"→ Updated HGNC ID from {transcript_hgnc_id} to {hgnc_id}",
        #                     file=sys.stderr,
        #                 )
        #                 tx["hgnc"] = hgnc_id
        #                 update_target["transcripts"].update({key: tx})
        #                 report.append(
        #                     ("transcript", "hgnc", key, transcript_hgnc_id, hgnc_id)
        #                 )
        # biotype_orig = list(sorted(tx.get("biotype", [])))
        # biotype = list(sorted(set(biotype_orig) | hgnc_to_biotype.get(_hgnc_id, set())))
        # if biotype_orig != biotype:
        #     tx["biotype"] = biotype
        #     update_target["transcripts"].update({key: tx})
        #     report.append(
        #         (
        #             "transcript",
        #             "biotype",
        #             key,
        #             ",".join(biotype_orig),
        #             ",".join(biotype),
        #         )
        #     )
        gene_symbol = tx.get("gene_symbol", None)
        if not gene_symbol:
            print("Transcript without gene symbol:", key, file=sys.stderr)
            if symbol := hgnc_to_symbol.get(_hgnc_id, None):
                symbol = symbol[0]
                print("→ Corresponding gene symbol:", symbol, file=sys.stderr)
                gene_symbol = symbol
                tx["gene_symbol"] = symbol
                update_target["transcripts"].update({key: tx})
            report.append(("transcript", "gene_symbol", key, None, gene_symbol))

        record = id_to_hgnc_record.get(_hgnc_id, {})
        if "Selenoproteins" in record.get("gene_group", []):
            print("Transcript is selenoprotein:", key, file=sys.stderr)
            biotype_orig = set(tx.get("biotype", []))
            biotype = list(sorted(biotype_orig | {"selenoprotein"}))
            tx["biotype"] = biotype
            update_target["transcripts"].update({key: tx})
            report.append(
                (
                    "transcript",
                    "biotype",
                    key,
                    ",".join(biotype_orig),
                    ",".join(biotype),
                )
            )

    return update_target, report


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    cdot, report = cdot_from_hgnc(
        snakemake.input.cdot,
        snakemake.input.hgnc,
    )

    with gzip.open(snakemake.output.cdot, "wt") as out:
        json.dump(cdot, out, indent=2)
        out.flush()
    df = pd.DataFrame(
        report, columns=["type", "what", "key", "hgnc_id_orig", "hgnc_id_new"]
    )
    df["kind"] = None
    df.loc[
        df["hgnc_id_orig"].notnull() & (df["hgnc_id_orig"] != df["hgnc_id_new"]), "kind"
    ] = "updated"
    df.loc[df["hgnc_id_orig"].isnull() & df["hgnc_id_new"].notnull(), "kind"] = (
        "fixed_missing"
    )
    df.loc[df["hgnc_id_orig"].isnull() & df["hgnc_id_new"].isnull(), "kind"] = (
        "still_missing"
    )
    df.to_csv(snakemake.output.report, sep="\t", index=False)
