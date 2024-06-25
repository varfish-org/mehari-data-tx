import sys
from io import BytesIO
from xml.etree.ElementTree import ElementTree

import requests

# url = r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmode=native&rettype=xml&id={transcript_ids}"
# ids = "&id=".join(snakemake.params.transcripts)
# ret = requests.get(url.format(transcript_ids=ids)).content
# with open(snakemake.output.xml, "wb") as f:
#     f.write(ret)
#
# f = BytesIO(ret)
with open(snakemake.output.xml, "r") as f:
    doc = ElementTree(file=f)

for e in doc.getroot().findall('./GBSeq'):
    intervals = e.findall(".//GBSeq_feature-table/GBFeature[GBFeature_key='exon']/GBFeature_intervals/GBInterval")
    for interval in intervals:
        accession = interval.findtext("GBInterval_accession")
        pos_from = interval.findtext("GBInterval_from")
        pos_to = interval.findtext("GBInterval_to")
        print(accession, pos_from, pos_to)
    cds = e.findall(".//GBFeature[GBFeature_key='CDS']/GBFeature_intervals/GBInterval")
    for interval in cds:
        cds_pos_from = interval.findtext("GBInterval_from")
        cds_pos_to = interval.findtext("GBInterval_to")
        print(f"{cds_pos_from=}..{cds_pos_to=}")

import gzip
with gzip.open(snakemake.output.cdot, "wt") as f:
    f.write("")
