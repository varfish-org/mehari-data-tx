#!/usr/bin/env python
"""Helper script to extract transcript identifier to tags from GRCh38.

The resulting file will be a TSV file and have the columns (1) tx_id,
(2) gene symbol, and (3) tags.  This can be used as the input for mehari
database building for applying the tags.
"""

import gzip
import json
import sys
import typing
from contextlib import redirect_stderr


def main(paths: typing.Iterable[str], file=sys.stdout):
    txs: dict[str, tuple[int, str, list[str]]] = {}
    for json_path in paths:
        print(f"processing {json_path}...", file=sys.stderr)

        print(f"- loading {json_path}", file=sys.stderr)
        if json_path.endswith(".gz"):
            with gzip.open(json_path, "rt") as inputf:
                json_data = json.load(inputf)
        else:
            with open(json_path, "rt") as inputf:
                json_data = json.load(inputf)

        print("- processing...", file=sys.stderr)
        for idx, tx in enumerate(json_data["transcripts"].values()):
            if idx > 0 and idx % 1000 == 0:
                print(f"  ... processed {idx} transcripts", file=sys.stderr)
            for genome_build in tx["genome_builds"].values():
                if genome_build.get("tag"):
                    tx_id, tx_version_ = tx["id"].split(".")
                    tx_version = int(tx_version_)
                    gene = tx["gene_name"]
                    tags = genome_build["tag"].split(",")
                    old_version, _, _ = txs.get(tx_id, (0, "", []))
                    if tx_version >= old_version:
                        txs[tx_id] = (tx_version, gene, tags)
        print(f"... done processing {json_path}", file=sys.stderr)
    print("finalizing...", file=sys.stderr)
    for tx_id, (tx_version, gene, tags) in sorted(txs.items()):
        lst = ",".join(sorted(tags))
        print(f"{tx_id}\t{gene}\t{lst}", file=file)
    print(file=file, flush=True)
    print("... done finalizing", file=sys.stderr)


with open(snakemake.output.mane_txs, "w") as out:
    with open(snakemake.log[0], "w") as log, redirect_stderr(log):
        main((snakemake.input.cdot,), file=out)
