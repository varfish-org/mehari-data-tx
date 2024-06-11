from contextlib import redirect_stderr

from snakemake import shell


def main():
    prefix = f"https://ftp.ncbi.nih.gov/refseq/{snakemake.params.species}/mRNA_Prot"

    # *.files.installed contains the list of files available
    url = f"{prefix}/{snakemake.params.species_name}.files.installed"

    # We're only interested in .rna.fna.gz files. For good measure, we sort them by their number.
    shell_cmd = (
        rf"""curl --silent {url} | grep ".rna.fna.gz" | cut -f2 | sort -t '.' -n -k2"""
    )

    # Download the respective file and append it to the output file.
    # Concatenating gzip files results in a valid (multi-)gzip file, so no need to do any (de/re-)compression.
    for part in shell(shell_cmd, iterable=True):
        shell(
            rf"""curl --silent '{prefix}/{part}' >> {snakemake.output.refseq_seq} 2>> {snakemake.log}"""
        )


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    main()
