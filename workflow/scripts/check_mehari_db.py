from contextlib import redirect_stderr

import polars as pl


def read_cdot_transcript_ids() -> set[str]:
    with open(snakemake.input.cdot_transcript_ids, "r") as f:
        cdot_transcript_ids = set(line.rstrip() for line in f)
        return cdot_transcript_ids


def read_cdot_hgnc_ids() -> set[str]:
    with open(snakemake.input.cdot_hgnc_ids, "r") as f:
        cdot_hgnc_ids = set(line.rstrip() for line in f)
        return cdot_hgnc_ids


def read_discarded_transcripts() -> pl.DataFrame:
    df = (
        pl.read_ndjson(
            snakemake.input.db_discarded,
            ignore_errors=True,
            schema={
                "type": pl.String,
                "value": pl.Struct(
                    {
                        "source": pl.String,
                        "reason": pl.String,
                        "id": pl.Struct(
                            {
                                "type": pl.String,
                                "value": pl.String,
                            }
                        ),
                        "gene_name": pl.String,
                        "tags": pl.String,
                    }
                ),
            },
        )
        .filter(pl.col("type") == "Discard")
        .unnest("value")
    )
    df = df.hstack(df.select(pl.col("id")).unnest("id").rename({"type": "value_type"}))
    return df


def get_discarded_transcript_ids(discarded: pl.DataFrame) -> set[str]:
    discarded_tx_ids = discarded.filter(pl.col("value_type").is_in(["TxId"])).select(
        pl.col("value")
    )
    return set(discarded_tx_ids.to_series().to_list())


def read_kept_transcript_ids() -> set[str]:
    with open(snakemake.input.kept_transcript_ids, "r") as f:
        return set(line.strip() for line in f)


def read_kept_gene_ids() -> set[str]:
    with open(snakemake.input.kept_hgnc_ids, "r") as f:
        return set(line.strip() for line in f)


def read_disease_genes() -> pl.DataFrame:
    df = pl.read_csv(snakemake.input.genes_to_disease, separator="\t")
    return df


def main():
    report = []

    cdot_tx_ids = read_cdot_transcript_ids()
    cdot_hgnc_ids = read_cdot_hgnc_ids()

    disease_gene_df = read_disease_genes()

    discarded = read_discarded_transcripts()
    stats: pl.DataFrame = discarded.group_by(["value_type", "reason"]).len()

    discarded_tx_ids = get_discarded_transcript_ids(discarded)

    kept_tx_ids = read_kept_transcript_ids()
    kept_hgnc_ids = read_kept_gene_ids()

    report.append(f"Number of transcripts in cdot:\t{len(cdot_tx_ids)}")
    report.append(f"Number of discarded transcripts:\t{len(discarded_tx_ids)}")
    report.append(f"Number of kept transcripts:\t{len(kept_tx_ids)}")

    num_kept = len(cdot_tx_ids & kept_tx_ids)
    num_discarded = len(cdot_tx_ids & discarded_tx_ids)

    report.append(f"Number of kept transcripts in cdot:\t{num_kept}")
    report.append(f"Number of discarded transcripts in cdot:\t{num_discarded}")

    valid = num_kept + num_discarded == len(cdot_tx_ids)

    for gene in disease_gene_df.iter_rows(named=True):
        hgnc_id = gene["hgnc_id"]
        if hgnc_id in kept_hgnc_ids:
            continue
        if hgnc_id not in cdot_hgnc_ids:
            report.append(f"Disease gene not in cdot:\t{hgnc_id}")
            continue
        hgnc_id = hgnc_id.lstrip("HGNC:")
        reason = discarded.filter(pl.col("value_type").is_in(["Hgnc"])).row(
            by_predicate=(pl.col("value") == hgnc_id),
            named=True,
        )["reason"]
        report.append(f"Discarded disease gene:\t{hgnc_id}\t{reason}")
        valid = False

    for transcript_id in cdot_tx_ids:
        tx_discarded = transcript_id in discarded_tx_ids
        tx_kept = transcript_id in kept_tx_ids
        if tx_discarded and tx_kept:
            report.append(f"Both kept and discarded:\t{transcript_id}")
            valid = False
        if not tx_discarded and not tx_kept:
            report.append(f"Neither kept nor discarded:\t{transcript_id}")
            valid = False

    discarded_mane = discarded.filter(
        pl.col("tags").str.contains("^.*Mane.*$"), pl.col("value_type").is_in(["Hgnc"])
    ).select(pl.col("value_type"), pl.col("value"), pl.col("gene_name"))
    for t, v, g in discarded_mane.iter_rows():
        report.append(f"Discarded MANE transcript:\t{t}:{v}\t{g}")
        valid = False

    NO_TRANSCRIPT_LEFT = ["NoTranscriptLeft"]
    investigate = discarded.filter(
        pl.col("value_type").is_in(["Hgnc"]),
        pl.col("reason").str.contains("NoTranscriptLeft"),
    ).select(
        pl.col("value_type"),
        pl.col("value"),
        pl.col("gene_name"),
        pl.col("reason")
        .str.split(" | ")
        .list.set_difference(NO_TRANSCRIPT_LEFT)
        .list.join(", "),
    )
    for t, v, g, r in investigate.rows():
        report.append(f"No Transcript left:\t{t}:{v}\t{g}\t{r}")
        # TODO check if this is a known issue, otherwise set valid = False

    report.append(f"Status:\t{'OK' if valid else 'ERROR'}")

    return stats, report


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    stats, report = main()
    stats.write_csv(snakemake.output.stats, separator="\t")
    with open(snakemake.output.report, "w") as f:
        f.write("\n".join(report))
