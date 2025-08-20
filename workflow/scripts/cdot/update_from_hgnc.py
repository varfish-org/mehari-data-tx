#!/usr/bin/env python
import gzip
import json
import sys
from collections import defaultdict
from contextlib import redirect_stderr
from copy import deepcopy
from dataclasses import dataclass
from typing import Type

import pandas as pd


@dataclass(frozen=True)
class HgncId:
    value: str


@dataclass(frozen=True)
class EntrezId:
    value: str


@dataclass(frozen=True)
class EnsemblGeneId:
    value: str


@dataclass(frozen=True)
class GeneSymbol:
    value: str


@dataclass(frozen=True)
class AliasSymbol:
    value: str


@dataclass(frozen=True)
class RefseqTranscriptId:
    value: str


@dataclass(frozen=True)
class EnsemblTranscriptId:
    value: str


@dataclass(frozen=True)
class GeneName:
    value: str


Identifier = EntrezId | EnsemblGeneId | GeneSymbol | AliasSymbol | RefseqTranscriptId | EnsemblTranscriptId | GeneName


def ids_with_optional_version(id_str: str | None, id_class: Type[Identifier]) -> list[Identifier]:
    """From a string identifier, creates dataclass instances for itself and its version-less form."""
    if not id_str:
        return []
    ids = [id_class(id_str)]
    if "." in id_str:
        ids.append(id_class(id_str.rsplit(".", 1)[0]))
    return ids


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


def build_identifier_hgnc_map(cdot: dict, hgnc: dict) -> tuple[defaultdict, defaultdict[HgncId, dict]]:
    """
    Builds a map of identifiers to HGNC Ids.
    Will use information from hgnc complete set first, then information from cdot json.
    """
    # identifier -> hgnc_id -> set of sources
    identifier_to_hgnc = defaultdict(lambda: defaultdict(set))
    hgnc_to_info = defaultdict(dict)  # store additional info such as symbol, biotype, etc., for each hgnc id

    # add an identifier to hgnc mapping, collecting all possible hgnc_ids and their sources
    def _add_to_map(identifier: Identifier, hgnc_id: HgncId, source: str):
        if identifier and identifier.value and hgnc_id:
            identifier_to_hgnc[identifier][hgnc_id].add(source)

    # process hgnc complete set
    hgnc_docs = hgnc.get("response", {}).get("docs", [])
    for record in hgnc_docs:
        hgnc_id = HgncId(record["hgnc_id"][5:])

        # store additional info
        hgnc_to_info[hgnc_id]["symbol"] = record.get("symbol")
        locus_type = record.get("locus_type")

        # update biotype based on locus_type
        if locus_type and LOCUS_TYPE_TO_BIOTYPE.get(locus_type) != "other":
            hgnc_to_info[hgnc_id].setdefault("biotype", set()).add(LOCUS_TYPE_TO_BIOTYPE[locus_type])
        if "Selenoproteins" in record.get("gene_group", []):
            hgnc_to_info[hgnc_id].setdefault("is_selenoprotein", True)

        if entrez_id := record.get("entrez_id"):
            _add_to_map(EntrezId(entrez_id), hgnc_id, "hgnc_entrez_id")
        if ensembl_gene_id := record.get("ensembl_gene_id"):
            for ident in ids_with_optional_version(ensembl_gene_id, EnsemblGeneId):
                _add_to_map(ident, hgnc_id, "hgnc_ensembl_gene_id")
        if symbol := record.get("symbol"):
            _add_to_map(GeneSymbol(symbol), hgnc_id, "hgnc_symbol")

        for accession in record.get("refseq_accession", []):
            for ident in ids_with_optional_version(accession, RefseqTranscriptId):
                _add_to_map(ident, hgnc_id, "hgnc_refseq_accession")
        for accession in record.get("mane_select", []):
            id_class = EnsemblTranscriptId if accession.startswith("ENS") else RefseqTranscriptId
            for ident in ids_with_optional_version(accession, id_class):
                _add_to_map(ident, hgnc_id, "hgnc_mane_select")
        for alias in record.get("alias_symbol", []):
            _add_to_map(AliasSymbol(alias), hgnc_id, "hgnc_alias_symbol")

    # process cdot json; the keys are the ncbi/ensembl gene ids
    for key, gene in cdot["genes"].items():
        if hgnc_id_str := gene.get("hgnc"):
            hgnc_id = HgncId(hgnc_id_str)
            id_class = EnsemblGeneId if key.startswith("ENSG") else EntrezId
            for ident in ids_with_optional_version(key, id_class):
                _add_to_map(ident, hgnc_id, "cdot_gene_key")
            if symbol := gene.get("gene_symbol"):
                _add_to_map(GeneSymbol(symbol), hgnc_id, "cdot_gene_symbol")

    # process cdot json; the keys are the ncbi/ensembl transcript ids
    for key, tx in cdot["transcripts"].items():
        if hgnc_id_str := tx.get("hgnc"):
            hgnc_id = HgncId(hgnc_id_str)
            id_class = EnsemblTranscriptId if key.startswith("ENST") else RefseqTranscriptId
            for ident in ids_with_optional_version(key, id_class):
                _add_to_map(ident, hgnc_id, "cdot_transcript_key")
            if gene_name := tx.get("gene_name"):
                _add_to_map(GeneName(gene_name), hgnc_id, "cdot_tx_gene_name")
            if gene_symbol := tx.get("gene_symbol"):
                _add_to_map(GeneSymbol(gene_symbol), hgnc_id, "cdot_tx_gene_symbol")
            if gene_version := tx.get("gene_version"):
                for ident in ids_with_optional_version(gene_version, EnsemblGeneId):
                    _add_to_map(ident, hgnc_id, "cdot_tx_gene_version")

    return identifier_to_hgnc, hgnc_to_info


def format_mapping_details(mappings: dict[HgncId, dict[Identifier, set[str]]]) -> str:
    parts = []
    for hgnc_id, sources in sorted(mappings.items(), key=lambda item: item[0].value):
        source_strs = []
        for ident, via in sorted(sources.items(), key=lambda item: str(item[0])):
            source_strs.append(f"{ident.__class__.__name__}('{ident.value}') (via {','.join(sorted(via))})")
        parts.append(f"HGNC:{hgnc_id.value} from {'; '.join(source_strs)}")
    return " | ".join(parts)


def update_cdot_from_map(cdot_to_update: dict, id_map: dict, info_map: dict[HgncId, dict]):
    """
    Go through cdot data and update it based on id_map and info_map.
    """
    report = []

    # update genes first
    for key, gene in cdot_to_update["genes"].items():
        hgnc_orig = gene.get("hgnc")
        identifiers: list[Identifier] = []
        id_class = EnsemblGeneId if key.startswith("ENSG") else EntrezId
        identifiers.extend(ids_with_optional_version(key, id_class))
        if symbol := gene.get("gene_symbol"):
            identifiers.append(GeneSymbol(symbol))
        if aliases := gene.get("aliases"):
            for alias in map(str.strip, aliases.split(",")):
                identifiers.append(AliasSymbol(alias))

        found_mappings: defaultdict[HgncId, defaultdict[Identifier, set[str]]] = defaultdict(lambda: defaultdict(set))
        for ident in filter(None, set(identifiers)):
            if ident in id_map:
                for hgnc_id, sources in id_map[ident].items():
                    found_mappings[hgnc_id][ident].update(sources)

        details = format_mapping_details(found_mappings) if found_mappings else ""
        hgnc_for_report, resolved_hgnc_id = hgnc_orig, None

        high_priority_sources = {"hgnc_ensembl_gene_id", "hgnc_entrez_id"}
        hgnc_by_direct_link = {
            hgnc_id
            for hgnc_id, sources in found_mappings.items()
            if any(src in high_priority_sources for src_set in sources.values() for src in src_set)
        }

        if len(hgnc_by_direct_link) == 1:
            single_id = hgnc_by_direct_link.pop()
            hgnc_for_report = single_id.value
            if hgnc_orig != single_id.value:
                gene["hgnc"] = single_id.value
            resolved_hgnc_id = single_id
        else:
            if len(found_mappings) == 1:
                single_id = next(iter(found_mappings.keys()))
                hgnc_for_report = single_id.value
                if hgnc_orig != single_id.value:
                    gene["hgnc"] = single_id.value
                resolved_hgnc_id = single_id
            elif len(found_mappings) > 1:
                conflicting_ids = {h.value for h in found_mappings.keys()}
                if hgnc_orig and hgnc_orig in conflicting_ids:
                    hgnc_for_report = hgnc_orig
                    resolved_hgnc_id = HgncId(hgnc_orig)
                    print(
                        f"CONFLICT for gene {key}: Original HGNC ID '{hgnc_orig}' is part of a conflict set {conflicting_ids}. Retaining original. Details: {details}",
                        file=sys.stderr,
                    )
                else:
                    conflicting_ids_str = ",".join(sorted(list(conflicting_ids)))
                    hgnc_for_report = conflicting_ids_str
                    print(f"CONFLICT for gene {key}: {details}", file=sys.stderr)

        report.append(("gene", "hgnc", key, hgnc_orig, hgnc_for_report, details))

        if resolved_hgnc_id:
            if not gene.get("gene_symbol") and (symbol := info_map.get(resolved_hgnc_id, {}).get("symbol")):
                gene["gene_symbol"] = symbol
                report.append(
                    ("gene", "gene_symbol", key, None, symbol, f"Added symbol for HGNC:{resolved_hgnc_id.value}")
                )

            biotype_orig = set(gene.get("biotype", []))
            added_biotypes = info_map.get(resolved_hgnc_id, {}).get("biotype", set())
            biotype_new_set = biotype_orig | added_biotypes
            if biotype_orig != biotype_new_set:
                biotype_new = sorted(list(biotype_new_set))
                gene["biotype"] = biotype_new
                report.append(
                    (
                        "gene",
                        "biotype",
                        key,
                        ",".join(sorted(list(biotype_orig))),
                        ",".join(biotype_new),
                        f"Added biotypes for HGNC:{resolved_hgnc_id.value}",
                    )
                )

    # update transcripts (similar logic to genes update above)
    for key, tx in cdot_to_update["transcripts"].items():
        hgnc_orig = tx.get("hgnc")
        identifiers: list[Identifier] = []
        id_class = EnsemblTranscriptId if key.startswith("ENST") else RefseqTranscriptId
        identifiers.extend(ids_with_optional_version(key, id_class))
        if tx_id := tx.get("id"):
            id_class = EnsemblTranscriptId if tx_id.startswith("ENST") else RefseqTranscriptId
            identifiers.extend(ids_with_optional_version(tx_id, id_class))
        if gene_name := tx.get("gene_name"):
            identifiers.append(GeneName(gene_name))
        if gene_symbol := tx.get("gene_symbol"):
            identifiers.append(GeneSymbol(gene_symbol))
        if gene_version := tx.get("gene_version"):
            identifiers.extend(ids_with_optional_version(gene_version, EnsemblGeneId))

        found_mappings: defaultdict[HgncId, defaultdict[Identifier, set[str]]] = defaultdict(lambda: defaultdict(set))
        for ident in filter(None, set(identifiers)):
            if ident in id_map:
                for hgnc_id, sources in id_map[ident].items():
                    found_mappings[hgnc_id][ident].update(sources)

        details = format_mapping_details(found_mappings) if found_mappings else ""
        hgnc_for_report, resolved_hgnc_id = hgnc_orig, None

        high_priority_tx_sources = {"hgnc_mane_select", "hgnc_refseq_accession"}
        hgnc_by_direct_tx_link = {
            hgnc_id
            for hgnc_id, sources in found_mappings.items()
            if any(src in high_priority_tx_sources for src_set in sources.values() for src in src_set)
        }

        if len(hgnc_by_direct_tx_link) == 1:
            single_id = hgnc_by_direct_tx_link.pop()
            hgnc_for_report = single_id.value
            if hgnc_orig != single_id.value:
                tx["hgnc"] = single_id.value
            resolved_hgnc_id = single_id
        else:
            if len(found_mappings) == 1:
                single_id = next(iter(found_mappings.keys()))
                hgnc_for_report = single_id.value
                if hgnc_orig != single_id.value:
                    tx["hgnc"] = single_id.value
                resolved_hgnc_id = single_id
            elif len(found_mappings) > 1:
                conflicting_ids = {h.value for h in found_mappings.keys()}
                if hgnc_orig and hgnc_orig in conflicting_ids:
                    hgnc_for_report = hgnc_orig
                    resolved_hgnc_id = HgncId(hgnc_orig)
                    print(
                        f"CONFLICT for transcript {key}: Original HGNC ID '{hgnc_orig}' is part of a conflict set {conflicting_ids}. Retaining original. Details: {details}",
                        file=sys.stderr,
                    )
                else:
                    conflicting_ids_str = ",".join(sorted(list(conflicting_ids)))
                    hgnc_for_report = conflicting_ids_str
                    print(f"CONFLICT for transcript {key}: {details}", file=sys.stderr)

        report.append(("transcript", "hgnc", key, hgnc_orig, hgnc_for_report, details))

        if resolved_hgnc_id:
            if not tx.get("gene_symbol") and (symbol := info_map.get(resolved_hgnc_id, {}).get("symbol")):
                tx["gene_symbol"] = symbol
                report.append(
                    ("transcript", "gene_symbol", key, None, symbol, f"Added symbol for HGNC:{resolved_hgnc_id.value}")
                )

            biotype_orig = set(tx.get("biotype", []))
            biotype_new_set = biotype_orig.copy()
            if info_map.get(resolved_hgnc_id, {}).get("is_selenoprotein"):
                biotype_new_set.add("selenoprotein")

            if biotype_orig != biotype_new_set:
                biotype_new = sorted(list(biotype_new_set))
                tx["biotype"] = biotype_new
                report.append(
                    (
                        "transcript",
                        "biotype",
                        key,
                        ",".join(sorted(list(biotype_orig))),
                        ",".join(biotype_new),
                        f"Added biotype for HGNC:{resolved_hgnc_id.value}",
                    )
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

    df = pd.DataFrame(report, columns=["type", "what", "key", "value_orig", "value_new", "details"])
    df["kind"] = "unchanged"

    df.loc[df["value_orig"].notnull() & (df["value_orig"] != df["value_new"]) & df["value_new"].notnull(), "kind"] = (
        "updated"
    )
    df.loc[df["value_orig"].isnull() & df["value_new"].notnull(), "kind"] = "fixed_missing"
    df.loc[df["value_orig"].notnull() & df["value_new"].isnull(), "kind"] = "removed"
    df.loc[df["value_orig"].isnull() & df["value_new"].isnull(), "kind"] = "still_missing"

    is_conflict = (df["what"] == "hgnc") & (df["value_new"].str.contains(",", na=False))
    df.loc[is_conflict, "kind"] = "conflict"

    is_retained_conflict = (
        (df["what"] == "hgnc")
        & (df["value_orig"] == df["value_new"])
        & (df["value_orig"].notnull())
        & (df["details"].str.contains(" | ", regex=False, na=False))
    )
    df.loc[is_retained_conflict, "kind"] = "conflict_retained_original"

    df = df[~((df["kind"] == "unchanged") & (df["what"] == "hgnc"))].copy()

    df.to_csv(snakemake.output.report, sep="\t", index=False)
