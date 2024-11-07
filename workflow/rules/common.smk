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
    params = config["sources"][alias]["cdot"]
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
    url = f"storage.googleapis.com/public-download-files/hgnc/archive/archive/{where}/json/hgnc_complete_set_{version}.json"
    return url


def cdot_input_mapping(wildcards: Wildcards) -> dict[str, str]:
    alias = get_alias(wildcards)
    if len(transcripts_to_fix_with_nuccore(wildcards)) > 0:
        cds = ".cds"
    else:
        cds = ""
    cdot_files = {
        "cdot": f"results/{alias}/cdot/{alias}.cdot.hgnc{cds}.json.gz",
    }
    if wildcards.source.lower() == "refseq":
        cdot_files.update(
            **{
                # "cdot_mt": f"results/{alias}/cdot/{alias}-from-ensembl.chrMT.json",
                "cdot_grafted": f"results/{alias}/cdot/{alias}.cdot.grafted.json.gz",
            }
        )
    return cdot_files


def ensembl_cdot(wildcards: Wildcards) -> str:
    alias_ensembl = f"{wildcards.assembly}-ensembl"
    return f"results/{alias_ensembl}/cdot/{alias_ensembl}.cdot.hgnc.json.gz"


def get_mehari_input(wildcards: Wildcards) -> dict[str, str]:
    alias = get_alias(wildcards)
    result = {
        "seqrepo_instance": f"results/{alias}/seqrepo/master",
        **cdot_input_mapping(wildcards),
    }
    if wildcards.assembly == "GRCh37":
        other = f"GRCh38-{wildcards.source}"
        result.update({"tags": f"results/{other}/cdot/{other}.tags.tsv"})
    return result


def get_mehari_cdot_param_string(wildcards: Wildcards, input: InputFiles) -> str:
    cdot_files = cdot_input_mapping(wildcards)
    expected_input_keys = input.keys()
    for key in cdot_files.keys():
        assert key in expected_input_keys
    return " ".join(f"--path-cdot-json {path}" for path in cdot_files.values())


def get_mehari_check_cdot_param_string(wildcards: Wildcards, input: InputFiles) -> str:
    cdot_files = cdot_input_mapping(wildcards)
    expected_input_keys = input.keys()
    for key in cdot_files.keys():
        assert key in expected_input_keys
    return " ".join(f"--cdot {path}" for path in cdot_files.values())


def transcripts_to_fix_with_nuccore(wildcards: Wildcards) -> list[str]:
    alias = get_alias(wildcards)
    return list(
        sorted(set(config["sources"][alias].get("fixes", {}).get("fix_cds", [])))
    )


def transcripts_to_lookup_ensembl_ids_for(wildcards: Wildcards) -> list[str]:
    alias = get_alias(wildcards)
    return list(
        sorted(
            set(config["sources"][alias].get("fixes", {}).get("add_from_ensembl", []))
        )
    )


def transcripts_to_fix_polyA(wildcards: Wildcards) -> list[str]:
    alias = get_alias(wildcards)
    return list(sorted(set(config["sources"][alias].get("fixes", {}).get("polyA", []))))


def parse_missing_accessions(_wildcards, input) -> list[str]:
    with open(input.missing_accessions) as file:
        accessions = {s.strip() for s in file}
    return list(sorted(accessions))


def get_genes_to_disease_download_url() -> str:
    url = "https://github.com/obophenotype/human-phenotype-ontology/releases/download/{release}/genes_to_disease.txt"
    release = config["human-phenotype-ontology"]["genes_to_disease"]["release"]
    return url.format(release=release)


def get_known_issues(wildcards: Wildcards) -> list[str]:
    alias = get_alias(wildcards)
    return config["sources"][alias].get("known_issues", [])
