from contextlib import redirect_stderr
import bioutils.seqfetcher


def main(missing_ids: str, file):
    with open(missing_ids, "r") as f:
        for line in f.readlines():
            accession = line.strip()
            seq = bioutils.seqfetcher.fetch_seq(accession)
            if seq:
                print(f">{accession}\n{seq}", file=file)


with open(snakemake.output.missing_fasta, "w") as out:
    with open(snakemake.log[0], "w") as log, redirect_stderr(log):
        main(snakemake.input.missing_txt, file=out)
