from typing import Callable

from snakemake.io import InputFiles, Wildcards


def get_ensembl_sequence_param(param: str) -> Callable[[Wildcards], str]:
    def inner(wildcards):
        return config["reference"][wildcards.alias]["ensembl"][param]

    return inner


def get_cdot_download_url(wildcards: Wildcards) -> str:
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


def get_hgnc_complete_set_download_url(_wildcards: Wildcards) -> str:
    version = config["hgnc"]["version"]
    if version[5:] in {"01-01", "04-01", "07-01", "10-01"}:
        where = "quarterly"
    else:
        where = "monthly"
    url = f"http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/{where}/json/hgnc_complete_set_{version}.json"
    return url


def _mehari_cdot_input(wildcards: Wildcards) -> dict[str, str]:
    alias = wildcards.alias
    cdot_files = {
        "cdot": f"results/transcripts/cdot/{alias}.hgnc.json.gz",
        "cdot_mt": "results/transcripts/cdot/GRCh38-ensembl.chrMT.json",
    }
    if alias == "GRCh38":
        cdot_files["cdot_graft"] = (
            "results/transcripts/cdot/GRCh38-ensembl.grafted.json.gz"
        )
    return cdot_files


def get_mehari_input(wildcards: Wildcards) -> dict[str, str]:
    alias = wildcards.alias
    seqrepo = wildcards.seqrepo
    result = {
        "seqrepo_instance": f"results/{seqrepo}/{alias}/master",
        **_mehari_cdot_input(wildcards),
    }
    if alias == "GRCh37":
        result.update({"mane_txs": "results/transcripts/cdot/GRCh37/mane-txs.tsv"})
    return result


def get_mehari_cdot_param_string(wildcards: Wildcards, input: InputFiles) -> str:
    cdot_files = _mehari_cdot_input(wildcards)
    expected_input_keys = input.keys()
    for key in cdot_files.keys():
        assert key in expected_input_keys
    return " ".join(f"--path-cdot-json {path}" for path in cdot_files.values())


def transcripts_to_fix_start_stop_codons_for(wildcards: Wildcards) -> set[str]:
    alias = wildcards.alias
    return set(config["transcripts"][alias]["fix_cds"])
