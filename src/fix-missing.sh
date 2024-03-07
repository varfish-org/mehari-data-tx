#!/usr/bin/env bash

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

set -x

# hotpatch broken seqrepo
sed -e 's/if aliases_cur.fetchone() is not None/if next(aliases_cur, None) is not None/' \
  /home/runner/micromamba-root/envs/mehari-data-tx/lib/python3.8/site-packages/biocommons/seqrepo/cli.py

mkdir -p $DATA_DIR/for-fix

{ grep "because it has no sequence" $DATA_DIR/pass-1/txs.bin.zst.report || true; } \
| cut -f 2 \
> $DATA_DIR/for-fix/missing.txt

wc -l $DATA_DIR/for-fix/missing.txt

xargs -a $DATA_DIR/for-fix/missing.txt \
seqrepo fetch-load -i master -n RefSeq \
>> $DATA_DIR/for-fix/missing.fasta
