import requests
import gzip

url = (
    r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
    r"db=nuccore"
    r"&retmode=native"
    r"&rettype=xml"
    r"&id={transcript_ids}"
)
ids = "&id=".join(snakemake.params.transcripts)
ret = requests.get(url.format(transcript_ids=ids)).content
with gzip.open(snakemake.output.xml, "wb") as f:
    f.write(ret)
