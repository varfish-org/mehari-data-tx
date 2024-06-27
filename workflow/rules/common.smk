from typing import Callable

from snakemake.io import InputFiles, Wildcards


def get_alias(wildcards: Wildcards) -> str:
    return f"{wildcards.assembly}-{wildcards.source}"


def genome_release(wildcards: Wildcards) -> str:
    return wildcards.assembly.lower()


def get_ensembl_sequence_param(param: str) -> Callable[[Wildcards], str]:
    def inner(wildcards):
        return config["reference"][wildcards.assembly]["ensembl"][param]

    return inner


def get_refseq_sequence_param(param: str) -> Callable[[Wildcards], str]:
    def inner(wildcards):
        return config["reference"][wildcards.assembly]["refseq"][param]

    return inner


def get_cdot_download_url(wildcards: Wildcards) -> str:
    alias = get_alias(wildcards)
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


def cdot_input_mapping(wildcards: Wildcards) -> dict[str, str]:
    alias = get_alias(wildcards)
    fix_order = config["fix-order"]
    cdot_files = {
        "cdot": f"results/{alias}/cdot/{alias}.cdot.fixed.hgnc.json.gz",
        "cdot_mt": f"results/{alias}/cdot/{alias}-from-ensembl.chrMT.json",
        "cdot_fixed": f"results/{alias}/cdot/{alias}.{fix_order}.json.gz",
    }
    return cdot_files


def ensembl_cdot(wildcards: Wildcards) -> str:
    alias_ensembl = wildcards.assembly + "-ensembl"
    return f"results/{alias_ensembl}/cdot/{alias_ensembl}.cdot.json.gz"


def get_mehari_input(wildcards: Wildcards) -> dict[str, str]:
    alias = get_alias(wildcards)
    seqrepo = wildcards.seqrepo
    result = {
        "seqrepo_instance": f"results/{alias}/{seqrepo}/master",
        **cdot_input_mapping(wildcards),
    }
    if alias == "GRCh37-refseq":
        result.update({"tags": f"results/{alias}/cdot/{alias}.tags.tsv"})
    return result


def get_mehari_cdot_param_string(wildcards: Wildcards, input: InputFiles) -> str:
    cdot_files = cdot_input_mapping(wildcards)
    expected_input_keys = input.keys()
    for key in cdot_files.keys():
        assert key in expected_input_keys
    return " ".join(f"--path-cdot-json {path}" for path in cdot_files.values())


def transcripts_to_fix_with_nuccore(wildcards: Wildcards) -> set[str]:
    alias = get_alias(wildcards)
    return set(config["transcripts"][alias]["fix_cds"])


def transcripts_to_lookup_ensembl_ids_for(wildcards: Wildcards) -> set[str]:
    alias = get_alias(wildcards)
    return set(config["transcripts"][alias]["add_from_ensembl"])


def missing_sequence_files(wildcards: Wildcards) -> list[str]:
    assembly = wildcards.assembly
    source = wildcards.source
    with checkpoints.detect_missing_sequences.get(
        assembly=assembly, source=source
    ).output["missing_txt"].open() as file:
        accessions = [s.strip() for s in file]
    return [
        f"results/{assembly}-{source}/mehari/seqrepo/missing/{accession}.fasta"
        for accession in accessions
    ]
