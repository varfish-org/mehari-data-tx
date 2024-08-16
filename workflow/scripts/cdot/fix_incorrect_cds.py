import gzip
import json
import logging
from collections import defaultdict, namedtuple
from contextlib import redirect_stderr
from datetime import datetime
from itertools import zip_longest
from typing import Any
from xml.etree.ElementTree import ElementTree

import pandas as pd
from snakemake.script import snakemake

Exon = namedtuple(
    "Exon",
    ["alt_start_i", "alt_end_i", "ord", "alt_cds_start_i", "alt_cds_end_i", "cigar"],
)


def update_cdot(cdot, doc):
    report = []
    exons = defaultdict(list[Exon])
    cds_positions = defaultdict(set)
    today = datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    logging.info(f"Getting information from NCBI nuccore ({today})")

    for e in doc.getroot().findall("./GBSeq"):
        intervals = e.findall(
            ".//GBSeq_feature-table/GBFeature[GBFeature_key='exon']/GBFeature_intervals/GBInterval"
        )

        for i, interval in enumerate(intervals):
            accession = interval.findtext("GBInterval_accession")
            pos_from = interval.findtext("GBInterval_from")
            pos_to = interval.findtext("GBInterval_to")
            exons[accession].append(Exon(-1, -1, i, int(pos_from), int(pos_to), None))

        cds = e.findall(
            ".//GBFeature[GBFeature_key='CDS']/GBFeature_intervals/GBInterval"
        )
        for interval in cds:
            accession = interval.findtext("GBInterval_accession")
            cds_pos_from = interval.findtext("GBInterval_from")
            cds_pos_to = interval.findtext("GBInterval_to")
            cds_positions[accession].add((int(cds_pos_from) - 1, int(cds_pos_to)))
        # TODO make use of
        #  <GBQualifier>
        #      <GBQualifier_name>codon_start</GBQualifier_name>
        #      <GBQualifier_value>${some_integer}</GBQualifier_value>
        #  </GBQualifier>
        #  to determine start codon position and adjust accordingly

    for accession in exons.keys() & cds_positions.keys():
        fix_exons(accession, cdot, exons, report)
        fix_cds(accession, cdot, cds_positions, report)

    return cdot, report


def fix_cds(
    accession: str,
    cdot: dict[str, Any],
    cds_positions: dict[str, set[tuple[int, int]]],
    report: list[tuple[str, str, Any, Any]],
):
    if (num := len(cds_positions[accession])) != 1:
        logging.warning(
            f"Expected exactly one CDS per transcript, got {num}. Likely a CDS join. "
            f"{accession=}:\n{cds_positions[accession]=}"
        )
    # assume start and stop codon are at from and to positions of the CDS
    if not cdot["transcripts"][accession].get("start_codon"):
        cdot["transcripts"][accession]["start_codon"] = None
    if not cdot["transcripts"][accession].get("stop_codon"):
        cdot["transcripts"][accession]["stop_codon"] = None
    report.append(
        (
            accession,
            "before",
            cdot["transcripts"][accession]["start_codon"],
            cdot["transcripts"][accession]["stop_codon"],
        )
    )
    cds_pos = next(iter(cds_positions[accession]))
    cdot["transcripts"][accession]["start_codon"] = cds_pos[0]
    cdot["transcripts"][accession]["stop_codon"] = cds_pos[1]
    report.append(
        (
            accession,
            "after",
            cdot["transcripts"][accession]["start_codon"],
            cdot["transcripts"][accession]["stop_codon"],
        )
    )


def fix_exons(
    accession: str,
    cdot: dict[str, Any],
    exons: list[Exon],
    report: list[tuple[str, str, Any, Any]],
):
    # Fix exon alignments as well if available; however, we have incomplete data here,
    # so the reference positions (alt_start_i, alt_end_i) are not available,
    # only the cds positions (alt_cds_start_i, alt_cds_end_i) are known.
    # The CIGAR string is also not available.
    if (tx := cdot["transcripts"].get(accession)) and (
        tx_exons := exons.get(accession)
    ):
        for genome_build in tx.get("genome_builds", []):
            current_exons = list(
                map(lambda e: Exon(*e), tx["genome_builds"][genome_build]["exons"])
            )

            def diffsum(a: list[Exon], b: list[Exon]):
                return sum(
                    aa.alt_cds_start_i != bb.alt_cds_start_i
                    or aa.alt_cds_end_i != bb.alt_cds_end_i
                    for aa, bb in zip_longest(
                        a, b, fillvalue=Exon(-1, -1, -1, -1, -1, None)
                    )
                )

            fwd_diff = diffsum(current_exons, tx_exons)
            rev_diff = diffsum(current_exons, tx_exons[::-1])
            if fwd_diff == 0 or rev_diff == 0:
                logging.info(
                    f"Exon cds positions for {accession} are the same as in current NCBI nuccore, skipping.",
                )
                continue

            # heuristic (because that information is not in the nuccore XML):
            # if the reverse order is more similar to the current alignment order, reverse the order
            if rev_diff < fwd_diff:
                tx_exons = [e for e in tx_exons[::-1]]
                for i, e in enumerate(tx_exons):
                    e = e._replace(ord=i)
                    tx_exons[i] = e

            report.append(
                (
                    accession,
                    "before",
                    tx["genome_builds"][genome_build]["exons"],
                    None,
                )
            )

            tx["genome_builds"][genome_build]["exons"] = list(map(list, tx_exons))

            report.append(
                (
                    accession,
                    "after",
                    tx["genome_builds"][genome_build]["exons"],
                    None,
                )
            )


def main(nuccore_xml_path: str, cdot_path: str) -> tuple[dict[str, Any], list]:
    with gzip.open(nuccore_xml_path, "r") as f:
        doc = ElementTree(file=f)

    with gzip.open(cdot_path, "r") as f:
        cdot = json.load(f)

    cdot, report = update_cdot(cdot, doc)
    return cdot, report


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    cdot, report = main(snakemake.input.xml, snakemake.input.cdot)

    with gzip.open(snakemake.output.cdot, "wt") as f:
        json.dump(cdot, f, indent=2)

    pd.DataFrame(
        report, columns=["accession", "change", "start_codon", "stop_codon"]
    ).to_csv(snakemake.output.report, sep="\t", index=False)
