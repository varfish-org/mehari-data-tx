import json
import os.path
import re
import sys
from contextlib import redirect_stderr
from io import BytesIO
from itertools import batched
from typing import Collection

import requests
from dinopy import FastaReader, FastaWriter
from snakemake import shell


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
    for batch in batched(accessions, max_ids):
        ids, versions = zip(*(accession.rsplit(".", 1) for accession in batch))
        data = {"ids": ids, "type": "cdna"}
        r = requests.post(server + ext, headers=headers, data=json.dumps(data))

        if not r.ok:
            raise ValueError(r)
        decoded = r.json()
        for accession, version in zip(ids, versions):
            if entry := decoded.get(accession):
                sequence = entry["seq"]
                if version != entry["version"]:
                    raise ValueError(
                        f"Version mismatch for {accession}: expected {version}, got {entry['version']}"
                    )
                yield sequence, f"{accession}.{version}"
            else:
                raise ValueError(f"Accession {accession} not found in response.")


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
        for entry in fetch_ncbi(ncbi_accessions):
            faw.write_entry(entry)
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
