import json
from contextlib import redirect_stderr

import polars as pl
import gzip


def cdot_transcript_ids() -> set[str]:
    with gzip.open(snakemake.input.cdot, "r") as f:
        cdot = json.load(f)
        cdot_transcript_ids = set(cdot["transcripts"].keys())
        return cdot_transcript_ids


def discarded_transcripts() -> pl.DataFrame:
    return (
        pl.scan_ndjson(snakemake.input.db_discarded, ignore_errors=True)
        .filter(pl.col("type") == "Discard")
        .unnest("value")
        .collect()
    )


def discarded_transcript_ids(discarded) -> set[str]:
    discarded_transcript_ids = (
        discarded.select(pl.col("id"))
        .unnest("id")
        .filter(pl.col("type").is_in(["TxId"]))
        .select(pl.col("value"))
    )
    return set(discarded_transcript_ids.to_series().to_list())


def kept_transcript_ids() -> set[str]:
    with open(snakemake.input.kept, "r") as f:
        return set(l.strip() for l in f)


def main():
    report = []

    cdot_tx_ids = cdot_transcript_ids()

    discarded = discarded_transcripts()
    stats: pl.DataFrame = (
        discarded.group_by(["source", "kind", "reason"])
        .len()
        .sort(["source", "kind", "len"])
    )

    discarded_tx_ids = discarded_transcript_ids(discarded)

    kept_tx_ids = kept_transcript_ids()

    report.append(f"Number of transcripts in cdot:\t{len(cdot_tx_ids)}")
    report.append(f"Number of discarded transcripts:\t{len(discarded_tx_ids)}")
    report.append(f"Number of kept transcripts:\t{len(kept_tx_ids)}")

    num_kept = len(cdot_tx_ids & kept_tx_ids)
    num_discarded = len(cdot_tx_ids & discarded_tx_ids)

    report.append(f"Number of kept transcripts in cdot:\t{num_kept}")
    report.append(f"Number of discarded transcripts in cdot:\t{num_discarded}")

    assert num_kept + num_discarded == len(cdot_tx_ids)

    for transcript_id in cdot_tx_ids:
        tx_discarded = transcript_id in discarded_tx_ids
        tx_kept = transcript_id in kept_tx_ids
        if tx_discarded and tx_kept:
            raise ValueError(f"Tx {transcript_id} is both kept and discarded")
        if not tx_discarded and not tx_kept:
            raise ValueError(f"Tx {transcript_id} is neither kept nor discarded")

    return stats, report


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    stats, report = main()
    stats.write_csv(snakemake.output.stats, separator="\t")
    with open(snakemake.output.report, "w") as f:
        f.write("\n".join(report))
