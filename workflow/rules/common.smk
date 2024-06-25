def get_ensembl_sequence_param(param):
    def inner(wildcards):
        return config["reference"][wildcards.alias]["ensembl"][param]

    return inner


def get_cdot_download_url(wildcards):
    alias = wildcards.alias
    params = config["transcripts"][alias]["cdot"]
    release = params["release"]
    custom = params.get("custom", None)
    source = params.get("source", "ensembl")
    if custom:
        url = f"https://github.com/SACGF/cdot/releases/download/data_v{release}/cdot-{release}.{custom}.json.gz"
    else:
        match source:
            case "ensembl":
                url = f"https://github.com/SACGF/cdot/releases/download/data_v{release}/cdot-{release}.{source}.{alias}.json.gz"
            case _:
                raise NotImplementedError(
                    f"Values other than 'ensembl' currently not supported ('{_}')"
                )
    return url


def get_hgnc_complete_set_download_url(_wildcards):
    version = config["hgnc"]["version"]
    if version[5:] in {"01-01", "04-01", "07-01", "10-01"}:
        where = "quarterly"
    else:
        where = "monthly"
    url = f"http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/{where}/json/hgnc_complete_set_{version}.json"
    return url


def get_mehari_input(wildcards):
    alias = wildcards.alias
    seqrepo = wildcards.seqrepo
    result = {
        "seqrepo_instance": f"results/{seqrepo}/{alias}/master",
        "cdot": f"results/transcripts/cdot/{alias}.json.gz",
        "cdot_hgnc": f"results/transcripts/cdot/{alias}.hgnc.json.gz",
        "cdot_mt": "results/transcripts/cdot/GRCh38-ensembl.chrMT.json",
    }
    if alias == "GRCh37":
        result.update({"mane_txs": "results/transcripts/cdot/GRCh37/mane-txs.tsv"})
    return result


def transcripts_to_fix_start_stop_codons_for(wildcards) -> set[str]:
    alias = wildcards.alias
    return set(config["transcripts"][alias]["fix_cds"])
