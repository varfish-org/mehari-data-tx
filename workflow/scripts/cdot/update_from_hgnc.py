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
    "RNA, long non-coding": "lncRNA",
    "RNA, micro": "other",
    "RNA, misc": "other",
    "RNA, ribosomal": "other",
    "RNA, small nuclear": "other",
    "RNA, small nucleolar": "other",
    "RNA, transfer": "other",
    "RNA, vault": "other",
    "RNA, Y": "Y_RNA",
    "T cell receptor gene": "other",
    "T cell receptor pseudogene": "other",
    "unknown": "other",
    "virus integration site": "other",
}


def build_identifier_hgnc_map(cdot: dict, hgnc: dict):
    """
    Builds a map of identifiers to HGNC Ids.
    Will use information from hgnc complete set first, then information from cdot json.
    """
    identifier_to_hgnc = defaultdict(str)  # mapping from any identifier to hgnc id
    hgnc_to_info = defaultdict(dict)  # store additional info such as symbol, biotype, etc., for each hgnc id

    # add an identifier to hgnc mapping, first write wins (TODO: add to set instead)
    def _add_to_map(identifier, hgnc_id):
        if identifier and hgnc_id and identifier not in identifier_to_hgnc:
            identifier_to_hgnc[identifier] = hgnc_id

    # process hgnc complete set
    hgnc_docs = hgnc.get("response", {}).get("docs", [])
    for record in hgnc_docs:
        hgnc_id = record["hgnc_id"][5:]

        # store additional info
        hgnc_to_info[hgnc_id]["symbol"] = record.get("symbol")
        locus_type = record.get("locus_type")

        # update biotype based on locus_type
        if locus_type and LOCUS_TYPE_TO_BIOTYPE.get(locus_type) != "other":
            hgnc_to_info[hgnc_id].setdefault("biotype", set()).add(LOCUS_TYPE_TO_BIOTYPE[locus_type])
        if "Selenoproteins" in record.get("gene_group", []):
            hgnc_to_info[hgnc_id].setdefault("is_selenoprotein", True)

        # map all known identifiers to this hgnc id
        # order matters: more reliable identifiers go first
        identifiers = [
            record.get("entrez_id"),
            record.get("ensembl_gene_id"),
            record.get("symbol"),
        ]
        identifiers.extend(record.get("refseq_accession", []))
        identifiers.extend(record.get("mane_select", []))
        identifiers.extend(record.get("alias_symbol", []))

        for identifier in identifiers:
            _add_to_map(identifier, hgnc_id)
            if isinstance(identifier, str) and "." in identifier:
                _add_to_map(identifier.rsplit(".", 1)[0], hgnc_id)

    # process cdot json; the keys are the ncbi/ensembl gene ids
    for key, gene in cdot["genes"].items():
        if hgnc_id := gene.get("hgnc"):
            _add_to_map(key, hgnc_id)
            if "." in key:
                _add_to_map(key.rsplit(".", 1)[0], hgnc_id)
            if symbol := gene.get("gene_symbol"):
                _add_to_map(symbol, hgnc_id)

    # process cdot json; the keys are the ncbi/ensembl transcript ids
    for key, tx in cdot["transcripts"].items():
        if hgnc_id := tx.get("hgnc"):
            _add_to_map(key, hgnc_id)
            if "." in key:
                _add_to_map(key.rsplit(".", 1)[0], hgnc_id)
            if gene_name := tx.get("gene_name"):
                _add_to_map(gene_name, hgnc_id)
            if gene_symbol := tx.get("gene_symbol"):
                _add_to_map(gene_symbol, hgnc_id)
            if gene_version := tx.get("gene_version"):
                _add_to_map(gene_version, hgnc_id)
                if "." in gene_version:
                    _add_to_map(gene_version.rsplit(".", 1)[0], hgnc_id)

    return identifier_to_hgnc, hgnc_to_info


def update_cdot_from_map(cdot_to_update: dict, id_map: dict, info_map: dict):
    """
    Go through cdot data and update it based on id_map and info_map.
    """
    report = []

    # update genes first
    for key, gene in cdot_to_update["genes"].items():
        hgnc_orig = gene.get("hgnc")

        # collect all known identifiers we want to look up in id_map
        identifiers = [key, gene.get("gene_symbol")]
        if "." in key:
            identifiers.append(key.rsplit(".", 1)[0])
        if aliases := gene.get("aliases"):
            identifiers.extend(map(str.strip, aliases.split(",")))

        # look up all identifiers in id_map
        found_ids = {id_map[i] for i in filter(None, set(identifiers)) if i in id_map}

        hgnc_for_report = hgnc_orig
        resolved_hgnc_id = None

        # check whether we have a unique hgnc id match
        if len(found_ids) == 1:
            # if so, modify cdot json accordingly
            single_id = found_ids.pop()
            hgnc_for_report = single_id
            if hgnc_orig != single_id:
                gene["hgnc"] = single_id
            resolved_hgnc_id = single_id
        elif len(found_ids) > 1:
            # otherwise: conflict detected. Do not modify cdot json and report conflict.
            conflicting_ids_str = ",".join(sorted(list(found_ids)))
            hgnc_for_report = conflicting_ids_str
            print(f"CONFLICT for gene {key}: Found multiple HGNC IDs: {hgnc_for_report}", file=sys.stderr)

        report.append(("gene", "hgnc", key, hgnc_orig, hgnc_for_report))

        # update symbol and biotype only if we have a resolved hgnc id
        if resolved_hgnc_id:
            if not gene.get("gene_symbol") and (symbol := info_map.get(resolved_hgnc_id, {}).get("symbol")):
                gene["gene_symbol"] = symbol
                report.append(("gene", "gene_symbol", key, None, symbol))

            biotype_orig = set(gene.get("biotype", []))
            added_biotypes = info_map.get(resolved_hgnc_id, {}).get("biotype", set())
            biotype_new_set = biotype_orig | added_biotypes
            if biotype_orig != biotype_new_set:
                biotype_new = sorted(list(biotype_new_set))
                gene["biotype"] = biotype_new
                report.append(("gene", "biotype", key, ",".join(sorted(list(biotype_orig))), ",".join(biotype_new)))

    # update transcripts (similar logic to genes update above)
    for key, tx in cdot_to_update["transcripts"].items():
        hgnc_orig = tx.get("hgnc")

        identifiers = [key, tx.get("id"), tx.get("gene_name"), tx.get("gene_symbol"), tx.get("gene_version")]
        if "." in key:
            identifiers.append(key.rsplit(".", 1)[0])
        if (gv := tx.get("gene_version")) and "." in gv:
            identifiers.append(gv.rsplit(".", 1)[0])

        found_ids = {id_map[i] for i in filter(None, set(identifiers)) if i in id_map}

        hgnc_for_report = hgnc_orig
        resolved_hgnc_id = None

        if len(found_ids) == 1:
            single_id = found_ids.pop()
            hgnc_for_report = single_id
            if hgnc_orig != single_id:
                tx["hgnc"] = single_id
            resolved_hgnc_id = single_id
        elif len(found_ids) > 1:
            conflicting_ids_str = ",".join(sorted(list(found_ids)))
            hgnc_for_report = conflicting_ids_str
            print(f"CONFLICT for transcript {key}: Found multiple HGNC IDs: {hgnc_for_report}", file=sys.stderr)

        report.append(("transcript", "hgnc", key, hgnc_orig, hgnc_for_report))

        if resolved_hgnc_id:
            if not tx.get("gene_symbol") and (symbol := info_map.get(resolved_hgnc_id, {}).get("symbol")):
                tx["gene_symbol"] = symbol
                report.append(("transcript", "gene_symbol", key, None, symbol))

            biotype_orig = set(tx.get("biotype", []))
            biotype_new_set = biotype_orig.copy()
            if info_map.get(resolved_hgnc_id, {}).get("is_selenoprotein"):
                biotype_new_set.add("selenoprotein")

            if biotype_orig != biotype_new_set:
                biotype_new = sorted(list(biotype_new_set))
                tx["biotype"] = biotype_new
                report.append(
                    ("transcript", "biotype", key, ",".join(sorted(list(biotype_orig))), ",".join(biotype_new))
                )

    return cdot_to_update, report


def cdot_from_hgnc(cdot_path: str, hgnc_path: str) -> tuple[dict, list[tuple[str, ...]]]:
    cdot_orig = load_json(cdot_path)
    hgnc_data = load_json(hgnc_path)
    id_map, info_map = build_identifier_hgnc_map(cdot_orig, hgnc_data)
    cdot_fixed = deepcopy(cdot_orig)
    cdot_updated, report = update_cdot_from_map(cdot_fixed, id_map, info_map)
    return cdot_updated, report


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    cdot, report = cdot_from_hgnc(
        snakemake.input.cdot,
        snakemake.input.hgnc,
    )

    with gzip.open(snakemake.output.cdot, "wt") as out:
        json.dump(cdot, out, indent=2)
        out.flush()
    df = pd.DataFrame(report, columns=["type", "what", "key", "hgnc_id_orig", "hgnc_id_new"])
    df["kind"] = "unchanged"

    df.loc[
        df["hgnc_id_orig"].notnull() & (df["hgnc_id_orig"] != df["hgnc_id_new"]) & df["hgnc_id_new"].notnull(), "kind"
    ] = "updated"
    df.loc[df["hgnc_id_orig"].isnull() & df["hgnc_id_new"].notnull(), "kind"] = "fixed_missing"
    df.loc[df["hgnc_id_orig"].notnull() & df["hgnc_id_new"].isnull(), "kind"] = "removed"
    df.loc[df["hgnc_id_orig"].isnull() & df["hgnc_id_new"].isnull(), "kind"] = "still_missing"

    is_conflict = df["hgnc_id_new"].str.contains(",", na=False)
    df.loc[is_conflict, "kind"] = "conflict"

    df.loc[df["kind"] == "unchanged", "hgnc_id_orig"] = df["hgnc_id_new"]
    df.to_csv(snakemake.output.report, sep="\t", index=False)
