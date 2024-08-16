import json
import os.path
import re
import sys
from contextlib import redirect_stderr
from io import BytesIO
from itertools import islice
from math import ceil
from time import sleep
from typing import Collection

import requests
from dinopy import FastaReader, FastaWriter
from snakemake.script import snakemake
from snakemake import shell


def batched(iterable, n):
    "Batch data into tuples of length n. The last batch may be shorter."
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch


def fetch_ncbi(accessions: Collection[str]):
    print(f"Fetching {len(accessions)} sequences from NCBI.", file=sys.stderr)
    accessions = list(sorted(accessions))
    dbfrom = "nuccore"
    db = "nucleotide"
    ids = ",".join(accessions)
    url_params = f"dbfrom={dbfrom}&db={db}&id={ids}"

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    epost_url = f"{base}epost.fcgi"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    response = requests.post(epost_url, data=url_params, headers=headers)
    output = response.text
    print(output, file=sys.stderr)

    if not output:
        raise ValueError("No output from eutils/esearch")

    # Parse WebEnv and QueryKey
    web = re.search(r"<WebEnv>(\S+)</WebEnv>", output).group(1)
    key = re.search(r"<QueryKey>(\d+)</QueryKey>", output).group(1)
    count = len(accessions)
    retmax = 500

    # Fetch sequences in batches of at most `retmax`
    for retstart in range(0, count, retmax):
        efetch_url = f"{base}efetch.fcgi?db=nucleotide&WebEnv={web}"
        efetch_url += f"&query_key={key}&retstart={retstart}"
        efetch_url += f"&retmax={retmax}&rettype=fasta&retmode=text"
        efetch_out = requests.get(efetch_url).text.strip().rstrip()
        if (
            "retstart" in efetch_out
            and "is+greater+than+number+of+records+available+in+history" in efetch_out
        ):
            break

        entries = efetch_out.replace("\n\n", "\n")
        # parse the output to make sure this really is FASTA data …
        far = FastaReader(BytesIO(entries.encode()))
        # … then yield the entries
        yield from far.entries()


def fetch_ensembl(accessions: Collection[str]):
    print(f"Fetching {len(accessions)} sequences from Ensembl.", file=sys.stderr)
    accessions = list(sorted(accessions))
    server = "https://rest.ensembl.org"
    ext = "/sequence/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    max_ids = 50
    batches = list(batched(accessions, max_ids))
    retry = False
    while len(batches) > 0:
        if not retry:
            batch = batches.pop(0)
        ids, versions = zip(*(accession.rsplit(".", 1) for accession in batch))
        data = {"ids": ids, "type": "cdna", "object_type": "transcript"}
        # To be on the safe side, we sleep for 1/15 seconds to have *at most* 15 requests per second.
        # For more thorough rate limiting behaviour, we should check the headers of the response.
        sleep(1.0 / 15.0)
        r = requests.post(server + ext, headers=headers, data=json.dumps(data))

        if not r.ok:
            if seconds := r.headers.get("Retry-After"):
                # n_requests = int(r.headers.get("X-RateLimit-Limit") or "0")
                # n_remaining = int(r.headers.get("X-RateLimit-Remaining") or "0")
                # period = int(r.headers.get("X-RateLimit-Period") or "0")
                # reset = int(r.headers.get("X-RateLimit-Reset") or "0")
                seconds = float(seconds)
                print(f"Waiting {seconds} seconds.", file=sys.stderr)
                retry = True
                sleep(int(ceil(seconds)))
                continue
            else:
                raise ValueError(r)

        retry = False
        decoded = r.json()
        if isinstance(decoded, list):
            decoded = {entry["id"]: entry for entry in decoded}
        for accession, version in zip(ids, versions):
            if entry := decoded.get(accession):
                sequence: str = entry["seq"]
                if int(version) != (remote_version := entry["version"]):
                    # TODO switch to biomart, filter on ensembl_transcript_id_version
                    # raise ValueError(
                    #    f"Version mismatch for {accession}: expected {version}, got {entry['version']}"
                    # )
                    print(
                        f"Version mismatch for {accession}: expected {version}, got {entry['version']}",
                        file=sys.stderr,
                    )
                    yield sequence.encode(), f"{accession}.{remote_version}".encode()
                yield sequence.encode(), f"{accession}.{version}".encode()
            else:
                print(f"Accession {accession} not found in response.", file=sys.stderr)
                # raise ValueError(f"Accession {accession} not found in response")


def main(accessions: Collection[str], fasta_path: str, fai_path: str):
    accessions = set(accessions)
    if os.path.exists(fasta_path) and os.path.exists(fai_path):
        with open(fai_path) as fai:
            for line in fai:
                accession = line.split("\t")[0]
                accessions.remove(accession)
    if not accessions:
        print("Nothing to do.", file=sys.stderr)
        if not os.path.exists(fasta_path):
            shell(f"touch {fasta_path}")
            shell(f"touch {fai_path}")
        exit(0)

    ensembl_accessions = {
        accession for accession in accessions if accession.startswith("ENS")
    }
    ncbi_accessions = accessions - ensembl_accessions

    with FastaWriter(fasta_path, append=True) as faw:
        if ncbi_accessions:
            for entry in fetch_ncbi(ncbi_accessions):
                faw.write_entry(entry)
        if ensembl_accessions:
            for entry in fetch_ensembl(ensembl_accessions):
                faw.write_entry(entry)

    # Create the index file, even if the FASTA file is empty
    if os.path.getsize(fasta_path) == 0:
        shell(f"touch {fai_path}")
    else:
        shell(f"samtools faidx {fasta_path}")


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    accessions = snakemake.params.missing_accessions
    fasta = snakemake.output.missing_fasta
    fai = snakemake.output.missing_fasta_fai
    main(accessions, fasta, fai)
