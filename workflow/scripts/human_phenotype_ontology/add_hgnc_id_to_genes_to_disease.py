import json
import sys
from contextlib import redirect_stderr

import pandas as pd


def main(genes_to_disease_path: str, hgnc_path: str) -> pd.DataFrame:
    with open(hgnc_path, "r") as f:
        hgnc = json.load(f)["response"]["docs"]
    hgnc = {
        f'NCBIGene:{gene["entrez_id"]}': gene["hgnc_id"]
        for gene in hgnc
        if "entrez_id" in gene and "hgnc_id" in gene
    }

    genes_to_disease = pd.read_csv(genes_to_disease_path, sep="\t")
    genes_to_disease["hgnc_id"] = genes_to_disease["ncbi_gene_id"].map(hgnc)

    discarded_genes = genes_to_disease[genes_to_disease["hgnc_id"].isnull()]
    print("Discarded genes without HGNC ID:", file=sys.stderr)
    print(discarded_genes["ncbi_gene_id"], file=sys.stderr)

    return genes_to_disease.dropna(subset=["hgnc_id"])


with redirect_stderr(open(snakemake.log[0], "w")):
    genes_to_disease = main(
        snakemake.input.genes_to_disease, snakemake.input.hgnc_complete_set
    )
    genes_to_disease.to_csv(snakemake.output.genes_to_disease, sep="\t", index=False)
