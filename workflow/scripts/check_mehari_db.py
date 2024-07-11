import sys
from contextlib import redirect_stderr

from polars import DataFrame
from snakemake.script import snakemake

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


def get_discarded_transcript_ids(discarded: DataFrame) -> set[str]:
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


def read_disease_genes() -> DataFrame:
    df = pl.read_csv(snakemake.input.genes_to_disease, separator="\t")
    return df


def check_known_issues(
    discarded: DataFrame, known_issues: set[str], report: list[str], valid: bool
):
    NO_TRANSCRIPT_LEFT: list[str] = ["NoTranscriptLeft"]
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
    for type_, type_id, gene_name, reason in investigate.rows():
        if type_id in known_issues:
            pass
        else:
            report.append(
                f"No Transcript left:\t{type_}:{type_id}\t{gene_name}\t{reason}"
            )
            valid = False
    return valid


def check_mane(discarded: DataFrame, report: list[str], valid: bool):
    discarded_mane = discarded.filter(
        pl.col("tags").str.contains("^.*Mane.*$"), pl.col("value_type").is_in(["Hgnc"])
    ).select(pl.col("value_type"), pl.col("value"), pl.col("gene_name"))
    for t, v, g in discarded_mane.iter_rows():
        report.append(f"Discarded MANE transcript:\t{t}:{v}\t{g}")
        valid = False
    return valid


def check_transcripts(
    cdot_tx_ids: set[str],
    discarded_tx_ids: set[str],
    kept_tx_ids: set[str],
    report: list[str],
    valid: bool,
):
    for transcript_id in cdot_tx_ids:
        tx_discarded = transcript_id in discarded_tx_ids
        tx_kept = transcript_id in kept_tx_ids
        if tx_discarded and tx_kept:
            report.append(f"Both kept and discarded:\t{transcript_id}")
            valid = False
        if not tx_discarded and not tx_kept:
            report.append(f"Neither kept nor discarded:\t{transcript_id}")
            valid = False
    return valid


def check_disease_genes(
    cdot_hgnc_ids: set[str],
    discarded: DataFrame,
    disease_gene_df: DataFrame,
    kept_hgnc_ids: set[str],
    report: list[str],
    valid: bool,
):
    for gene in disease_gene_df.iter_rows(named=True):
        hgnc_id = gene["hgnc_id"]
        if hgnc_id in kept_hgnc_ids:
            continue

        hgnc_id = hgnc_id.lstrip("HGNC:")
        if hgnc_id not in cdot_hgnc_ids:
            report.append(
                f"Disease gene not in cdot:\tHGNC:{hgnc_id}\t{gene['gene_symbol']}"
            )
            continue

        reason = discarded.filter(pl.col("value_type").is_in(["Hgnc"])).row(
            by_predicate=(pl.col("value") == hgnc_id),
            named=True,
        )["reason"]
        if "NoTranscripts" in reason:
            print(
                f"Discarded disease gene  :\tHGNC:{hgnc_id}\t{gene['gene_symbol']} (No Transcripts to begin with)",
                file=sys.stderr,
            )
        else:
            report.append(
                f"Discarded disease gene  :\tHGNC:{hgnc_id}\t{gene['gene_symbol']}\t{reason}"
            )
            valid = False
    return valid


def main(known_issues: set[str]):
    report: list[str] = []

    cdot_tx_ids: set[str] = read_cdot_transcript_ids()
    cdot_hgnc_ids: set[str] = read_cdot_hgnc_ids()

    disease_gene_df: DataFrame = read_disease_genes()

    discarded: DataFrame = read_discarded_transcripts()
    stats: DataFrame = discarded.group_by(["value_type", "reason"]).len()

    discarded_tx_ids: set[str] = get_discarded_transcript_ids(discarded)

    kept_tx_ids: set[str] = read_kept_transcript_ids()
    kept_hgnc_ids: set[str] = read_kept_gene_ids()

    report.append(f"Number of transcripts in cdot:\t{len(cdot_tx_ids)}")
    report.append(f"Number of discarded transcripts:\t{len(discarded_tx_ids)}")
    report.append(f"Number of kept transcripts:\t{len(kept_tx_ids)}")

    num_kept = len(cdot_tx_ids & kept_tx_ids)
    num_discarded = len(cdot_tx_ids & discarded_tx_ids)

    report.append(f"Number of kept transcripts in cdot:\t{num_kept}")
    report.append(f"Number of discarded transcripts in cdot:\t{num_discarded}")

    valid: bool = num_kept + num_discarded == len(cdot_tx_ids)

    valid = check_disease_genes(
        cdot_hgnc_ids, discarded, disease_gene_df, kept_hgnc_ids, report, valid
    )
    valid = check_transcripts(cdot_tx_ids, discarded_tx_ids, kept_tx_ids, report, valid)
    valid = check_mane(discarded, report, valid)
    valid = check_known_issues(discarded, known_issues, report, valid)

    report.append(f"Status:\t{'OK' if valid else 'ERROR'}")

    return stats, report


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    known_issues: set[str] = set(snakemake.params.known_issues)
    stats, report = main(known_issues=known_issues)
    stats.write_csv(snakemake.output.stats, separator="\t")
    with open(snakemake.output.report, "w") as f:
        f.write("\n".join(report))
