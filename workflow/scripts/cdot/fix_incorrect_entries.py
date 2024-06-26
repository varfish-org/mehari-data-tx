import json
import logging
from collections import defaultdict
from contextlib import redirect_stderr
from typing import Any
from xml.etree.ElementTree import ElementTree
import gzip

import pandas as pd


def update_cdot(cdot, doc):
    report = []
    exons = defaultdict(list)
    cds_positions = defaultdict(set)
    for e in doc.getroot().findall("./GBSeq"):
        intervals = e.findall(
            ".//GBSeq_feature-table/GBFeature[GBFeature_key='exon']/GBFeature_intervals/GBInterval"
        )

        for interval in intervals:
            accession = interval.findtext("GBInterval_accession")
            pos_from = interval.findtext("GBInterval_from")
            pos_to = interval.findtext("GBInterval_to")
            exons[accession].append((int(pos_from) - 1, int(pos_to)))

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
        # skip fixing exons for now until we know what cdot needs
        # cdot["transcripts"][accession]["genome_builds"][genome_build]["exons"]

        if (num := len(cds_positions[accession])) != 1:
            logging.warning(
                f"Expected exactly one CDS per transcript, got {num}. Likely a CDS join. "
                f"{accession=}:\n{cds_positions[accession]=}"
            )
        # assume start and stop codon are at from and to positions of the CDS
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
    return cdot, report


def main(
    genome_build: str, nuccore_xml_path: str, cdot_path: str
) -> tuple[dict[str, Any], list]:
    with gzip.open(nuccore_xml_path, "r") as f:
        doc = ElementTree(file=f)

    with gzip.open(cdot_path, "r") as f:
        cdot = json.load(f)

    if genome_build == "GRCh38":
        cdot, report = update_cdot(cdot, doc)
        return cdot, report
    return cdot, []


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    cdot, report = main(
        snakemake.wildcards.alias, snakemake.input.xml, snakemake.input.cdot
    )

    with gzip.open(snakemake.output.cdot, "wt") as f:
        json.dump(cdot, f, indent=2)

    pd.DataFrame(
        report, columns=["accession", "change", "start_codon", "stop_codon"]
    ).to_csv(snakemake.output.report, sep="\t", index=False)
