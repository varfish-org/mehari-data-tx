#!/usr/bin/env python
"""Set the MANE tags of a GRCh38 RefSeq cdot JSON from the MANE summary.

cdot lacks some MANE tags, e.g. on non-coding MANE transcripts and on MANE transcripts
whose sequence differs from the genome. The summary lists all MANE transcripts, so its
MANE status replaces the MANE tags of cdot. Transcripts match by accession and version.
"""

import csv
import gzip
import json
import sys
from contextlib import redirect_stderr

MANE_TAGS = ("MANE Select", "MANE Plus Clinical")


def check_versions(path: str, refseq_release: str):
    """Fail unless the MANE release belongs to the RefSeq release."""
    with open(path) as file:
        versions = dict(line.rstrip("\n").split("\t", 1) for line in file if "\t" in line)
    annotation = versions.get("NCBI RefSeq Annotation Release", "")
    if not annotation.endswith(f"-{refseq_release}"):
        raise ValueError(
            f"MANE {versions.get('MANE Version')} belongs to {annotation}, not to {refseq_release}"
        )


def main(cdot: dict, summary_path: str) -> dict:
    with gzip.open(summary_path, "rt") as file:
        mane = {row["RefSeq_nuc"]: row for row in csv.DictReader(file, delimiter="\t")}
    for tx_id, tx in cdot["transcripts"].items():
        alignment = tx["genome_builds"].get("GRCh38")
        if alignment is None:
            continue
        tags = [tag for tag in (alignment.get("tag") or "").split(",") if tag]
        new_tags = [tag for tag in tags if tag not in MANE_TAGS]
        if tx_id in mane:
            new_tags.append(mane[tx_id]["MANE_status"])
        if sorted(new_tags) != sorted(tags):
            print(f"Updated transcript {tx_id} tags from {tags} to {new_tags}", file=sys.stderr)
            if new_tags:
                alignment["tag"] = ",".join(new_tags)
            else:
                del alignment["tag"]
    for tx_id in sorted(mane.keys() - cdot["transcripts"].keys()):
        row = mane[tx_id]
        print(
            f"Not in the annotation: {tx_id} ({row['MANE_status']} of {row['symbol']} on {row['GRCh38_chr']})",
            file=sys.stderr,
        )
    return cdot


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    check_versions(snakemake.input.versions, snakemake.params.refseq_release)
    with gzip.open(snakemake.input.cdot, "rt") as file:
        cdot = json.load(file)
    cdot = main(cdot, snakemake.input.summary)
    with gzip.open(snakemake.output.cdot, "wt") as out:
        json.dump(cdot, out, indent=2)
