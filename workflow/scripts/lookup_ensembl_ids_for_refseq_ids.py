from contextlib import redirect_stderr

import requests


def main(refseq_ids: list[str]):
    refseq_ids = ",".join(s.rsplit(".", 1)[0] for s in refseq_ids)

    xml = rf"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
        <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
            <Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y"/>
            <Filter name = "refseq_mrna" value = "{refseq_ids}"/>\\
            <Attribute name = "refseq_mrna" />
            <Attribute name = "ensembl_gene_id" />\
            <Attribute name = "ensembl_gene_id_version" />\
            <Attribute name = "ensembl_transcript_id" />
            <Attribute name = "ensembl_transcript_id_version" />
            <Attribute name = "hgnc_symbol" />
            <Attribute name = "hgnc_id" />
        </Dataset>
    </Query>"""

    biomart_url = rf"http://www.ensembl.org/biomart/martservice?query={xml}"
    response = requests.get(biomart_url)
    return response


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    response = main(snakemake.params.refseq_ids)
    with open(snakemake.output.tsv, "wt") as f:
        f.write(response.content.decode())
