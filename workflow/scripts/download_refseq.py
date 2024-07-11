from contextlib import redirect_stderr
from snakemake.script import snakemake

from snakemake import shell


def main(installed_files: list[str]):
    prefix = f"https://ftp.ncbi.nih.gov/refseq/{snakemake.params.species}/mRNA_Prot"

    # Download the respective file and append it to the output file.
    # Concatenating gzip files results in a valid (multi-)gzip file, so no need to do any (de/re-)compression.
    for part in installed_files:
        shell(
            rf"""curl --silent '{prefix}/{part}' >> {snakemake.output.refseq_seq} 2>> {snakemake.log}"""
        )


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    with open(snakemake.input.refseq_installed, "rt") as file:
        lines = [line.strip() for line in file]
        files = [line.split() for line in lines if line.endswith(".rna.fna.gz")]
        files = sorted(files, key=lambda x: int(x[1].split(".")[1]))
        _checksums, files = zip(*files)
    main(files)
