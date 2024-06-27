from contextlib import redirect_stderr
from datetime import timedelta

import bioutils.seqfetcher
from ratelimit import limits, sleep_and_retry


@sleep_and_retry
@limits(calls=10, period=timedelta(seconds=1).total_seconds())
def get_fasta(accession: str) -> str | None:
    try:
        seq = bioutils.seqfetcher.fetch_seq(accession)
        if seq:
            return f">{accession}\n{seq}"
    except Exception as e:
        print(f"Failed to fetch {accession}: {e}", file=sys.stderr)
    return None


def main(missing_ids: str, file):
    with open(missing_ids, "r") as f:
        for line in f.readlines():
            accession = line.strip()
            fasta = get_fasta(accession)
            if fasta:
                print(fasta, file=file)


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    with open(snakemake.output.missing_fasta, "w") as out:
        main(snakemake.input.missing_txt, file=out)
