from contextlib import redirect_stderr

import requests
import gzip


def main(transcript_ids: list[str]) -> bytes:
    url = (
        r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        r"db=nuccore"
        r"&retmode=native"
        r"&rettype=xml"
        r"&id={transcript_ids}"
    )
    ids = "&id=".join(transcript_ids)
    ret = requests.get(url.format(transcript_ids=ids)).content
    return ret


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    if len(snakemake.params.transcripts) > 0:
        ret = main(snakemake.params.transcripts)
    else:
        ret = b"<GBSet></GBSet>"
    with gzip.open(snakemake.output.xml, "wb") as f:
        f.write(ret)
