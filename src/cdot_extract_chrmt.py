#!/usr/bin/env python
"""Helper script to alignments to chrMT from CDOT JSON.

This is used to graft ENSEMBL chrMT transcripts to the RefSeq transcripts
database.
"""

import gzip
import json
import sys
import typing

#: Contig name for chrMT
CHRMT_CONTIG = "NC_012920.1"


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <json_path>", file=sys.stderr)
        return 1
    json_path = sys.argv[1]

    #: Transcript IDs (with version) to keep.
    tx_ids: typing.Set[str] = set()
    #: Gene IDs (with version) tok eep.
    gene_ids: typing.Set[str] = set()

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
            if genome_build["contig"] == CHRMT_CONTIG:
                # print(json.dump(tx, sys.stderr, indent=2))
                # print("\n", file=sys.stderr)
                tx_ids.add(tx["id"])
                gene_ids.add(tx["gene_version"])
    print(f"... done processing {json_path}", file=sys.stderr)

    print("reducing to chrMT...", file=sys.stderr)
    for gene_id in list(json_data["genes"].keys()):
        if gene_id not in gene_ids:
            json_data["genes"].pop(gene_id)
    for tx_id in list(json_data["transcripts"].keys()):
        if tx_id not in tx_ids:
            json_data["transcripts"].pop(tx_id)
    print("... done reducing to chrMT", file=sys.stderr)

    json.dump(json_data, sys.stdout, indent=2)
    print(file=sys.stdout)


if __name__ == "__main__":
    sys.exit(main())
