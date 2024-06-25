import json
from collections import defaultdict
from xml.etree.ElementTree import ElementTree
import gzip

genome_build = snakemake.wildcards.alias

with gzip.open(snakemake.input.xml, "r") as f:
    doc = ElementTree(file=f)

with gzip.open(snakemake.input.cdot, "r") as f:
    cdot = json.load(f)


def update_cdot(cdot, doc):
    exons = defaultdict(list)
    cds_positions = defaultdict(list)
    for e in doc.getroot().findall('./GBSeq'):
        intervals = e.findall(".//GBSeq_feature-table/GBFeature[GBFeature_key='exon']/GBFeature_intervals/GBInterval")

        for interval in intervals:
            accession = interval.findtext("GBInterval_accession")
            pos_from = interval.findtext("GBInterval_from")
            pos_to = interval.findtext("GBInterval_to")
            exons[accession].append((int(pos_from), int(pos_to)))

        cds = e.findall(".//GBFeature[GBFeature_key='CDS']/GBFeature_intervals/GBInterval")
        for interval in cds:
            accession = interval.findtext("GBInterval_accession")
            cds_pos_from = interval.findtext("GBInterval_from")
            cds_pos_to = interval.findtext("GBInterval_to")
            cds_positions[accession].append((int(cds_pos_from), int(cds_pos_to)))

    for accession in exons.keys() & cds_positions.keys():
        # skip fixing exons for now until we know what cdot needs
        # cdot["transcripts"][accession]["genome_builds"][genome_build]["exons"]

        assert len(cds_positions[accession]) == 1, "Expected exactly one CDS per transcript, got more"
        # assume start and stop codon are at from and to positions of the CDS
        cdot["transcripts"][accession]["start_codon"] = cds_positions[accession][0][0]
        cdot["transcripts"][accession]["stop_codon"] = cds_positions[accession][0][1]
    return cdot


if genome_build == "GRCh38":
    cdot = update_cdot(cdot, doc)

with gzip.open(snakemake.output.cdot, "wt") as f:
    json.dump(cdot, f, indent=2)
