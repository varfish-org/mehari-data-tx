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


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    with open(snakemake.output.missing_fasta, "w") as out:
        fasta = get_fasta(snakemake.params.accession)
        if fasta:
            print(fasta, file=file)
