from contextlib import redirect_stderr

import polars as pl


def cdot_transcript_ids() -> set[str]:
    with open(snakemake.input.cdot_tx_ids, "r") as f:
        cdot_transcript_ids = set(line.rstrip() for line in f)
        return cdot_transcript_ids


def discarded_transcripts() -> pl.DataFrame:
    return (
        pl.scan_ndjson(snakemake.input.db_discarded, ignore_errors=True)
        .filter(pl.col("type") == "Discard")
        .unnest("value")
        .collect()
    )


def discarded_transcript_ids(discarded: pl.DataFrame) -> set[str]:
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
    stats: pl.DataFrame = discarded.group_by(["source", "kind", "reason"]).len()

    discarded_tx_ids = discarded_transcript_ids(discarded)

    kept_tx_ids = kept_transcript_ids()

    report.append(f"Number of transcripts in cdot:\t{len(cdot_tx_ids)}")
    report.append(f"Number of discarded transcripts:\t{len(discarded_tx_ids)}")
    report.append(f"Number of kept transcripts:\t{len(kept_tx_ids)}")

    num_kept = len(cdot_tx_ids & kept_tx_ids)
    num_discarded = len(cdot_tx_ids & discarded_tx_ids)

    report.append(f"Number of kept transcripts in cdot:\t{num_kept}")
    report.append(f"Number of discarded transcripts in cdot:\t{num_discarded}")

    valid = num_kept + num_discarded == len(cdot_tx_ids)
    for transcript_id in cdot_tx_ids:
        tx_discarded = transcript_id in discarded_tx_ids
        tx_kept = transcript_id in kept_tx_ids
        if tx_discarded and tx_kept:
            report.append(f"Both kept and discarded:\t{transcript_id}")
            valid = False
        if not tx_discarded and not tx_kept:
            report.append(f"Neither kept nor discarded:\t{transcript_id}")
            valid = False
    report.append(f"Status:\t{'OK' if valid else 'ERROR'}")

    return stats, report


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    stats, report = main()
    stats.write_csv(snakemake.output.stats, separator="\t")
    with open(snakemake.output.report, "w") as f:
        f.write("\n".join(report))
