from contextlib import redirect_stderr
from io import StringIO

import pandas as pd
import requests
from snakemake.script import snakemake

HEADER = {
    "refseq_mrna": "RefSeq mRNA ID",
    "ensembl_gene_id": "Gene stable ID",
    "ensembl_gene_id_version": "Gene stable ID version",
    "ensembl_transcript_id": "Transcript stable ID",
    "ensembl_transcript_id_version": "Transcript stable ID version",
    "hgnc_id": "HGNC ID",
}

CONTIGS = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "MT",
    "X",
    "Y",
]


def main(refseq_ids: set[str]):
    refseq_ids = ",".join(s.rsplit(".", 1)[0] for s in refseq_ids)
    attributes = "".join(f'<Attribute name = "{name}" />' for name in HEADER.keys())
    contigs = ",".join(CONTIGS)
    xml = rf"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
        <Dataset name="hsapiens_gene_ensembl" interface="default" >
            <Filter name="chromosome_name" value="{contigs}"/>
            <Filter name="refseq_mrna" value="{refseq_ids}"/>\\
            {attributes}
        </Dataset>
    </Query>"""

    biomart_url = rf"http://www.ensembl.org/biomart/martservice?query={xml}"
    response = requests.get(biomart_url)
    return response


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    refseq_ids = set(snakemake.params.additional_accessions)
    with open(snakemake.input.accessions, "rt") as f:
        for line in f:
            refseq_ids.add(line.strip())
    columns = list(HEADER.values())
    if snakemake.wildcards.source.lower() == "refseq":
        response = main(refseq_ids).content.decode()
        try:
            df = pd.read_csv(StringIO(response), sep="\t", header=0, usecols=columns)
        except Exception as e:
            raise ValueError(
                f"Expected TSV response from biomart, got:\n{response}"
            ) from e
    else:
        # Built an empty file for the case of Ensembl / non-refseq
        response = "\t".join(columns) + "\n"
        df = pd.read_csv(StringIO(response), sep="\t", header=0, usecols=columns)

    df.to_csv(snakemake.output.tsv, sep="\t", index=False)
