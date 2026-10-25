from typing import Callable

from snakemake.iocontainers import InputFiles, Wildcards


def get_alias(wildcards: Wildcards) -> str:
    return f"{wildcards.assembly}-{wildcards.source}"


def source_kind(wildcards: Wildcards) -> str:
    """The source without its release: 'ensembl' or 'refseq'."""
    return wildcards.source.split("-", 1)[0]


def source_release(wildcards: Wildcards) -> str:
    """The Ensembl or RefSeq annotation release of the source."""
    return wildcards.source.split("-", 1)[1]


def get_source_config_key(wildcards: Wildcards) -> str:
    """Key in config["sources"]; all releases of a source share one entry."""
    return f"{wildcards.assembly}-{source_kind(wildcards)}"


def releases(assembly: str, kind: str) -> list:
    return config["reference"][assembly][kind]["releases"]


def newest_release(assembly: str, kind: str):
    """The highest Ensembl release, or the last listed RefSeq release."""
    if kind == "ensembl":
        return max(releases(assembly, kind))
    return releases(assembly, kind)[-1]


def newest_alias(assembly: str, kind: str) -> str:
    return f"{assembly}-{kind}-{newest_release(assembly, kind)}"


def used_alias(wildcards: Wildcards, kind: str) -> str:
    """The alias of the given kind whose sequences and annotation a source uses.

    That is the source itself, or the newest release of the other source.
    """
    if source_kind(wildcards) == kind:
        return get_alias(wildcards)
    return newest_alias(wildcards.assembly, kind)


def reference_path(kind: str, suffix: str) -> Callable[[Wildcards], str]:
    def inner(wildcards):
        alias = used_alias(wildcards, kind)
        return f"results/{alias}/reference/{alias}.{suffix}"

    return inner


def all_aliases() -> list[str]:
    """Aliases of all databases built from a single source."""
    return [
        f"{assembly}-{kind}-{release}"
        for assembly in config["reference"]
        for kind in ("ensembl", "refseq")
        for release in releases(assembly, kind)
    ]


def merged_alias(assembly: str) -> str:
    return (
        f"{newest_alias(assembly, 'ensembl')}"
        f"-and-refseq-{newest_release(assembly, 'refseq')}"
    )


def genome_assembly(wildcards: Wildcards) -> str:
    return wildcards.assembly.lower()


def genome_assembly_version_parameter(_wildcards: Wildcards) -> str:
    # TODO determine e.g. patch version for GRCh37/38
    return ""


def transcript_source(wildcards: Wildcards) -> str:
    return source_kind(wildcards)


def transcript_source_version_parameter(wildcards: Wildcards) -> str:
    return f"--transcript-source-version {source_release(wildcards)}"


def cdot_version(wildcards: Wildcards) -> str:
    return config["sources"][get_source_config_key(wildcards)]["cdot"]["release"]


def get_ensembl_sequence_param(param: str) -> Callable[[Wildcards], str]:
    def inner(wildcards):
        return config["reference"][wildcards.assembly]["ensembl"][param]

    return inner


def get_refseq_sequence_url(wildcards: Wildcards) -> str:
    url = config["reference"][wildcards.assembly]["refseq"]["fasta"]
    return url.format(refseq_release=wildcards.refseq_release)


def get_cdot_download_url(wildcards: Wildcards) -> str:
    alias = get_alias(wildcards)
    params = config["sources"][get_source_config_key(wildcards)]["cdot"]
    release = params["release"]
    custom = params.get("custom", None)
    source = params.get("source", "ensembl")
    if custom:
        custom = custom.format(
            **{f"{source_kind(wildcards)}_release": source_release(wildcards)}
        )
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
    archive = config["hgnc"]["archive"]
    url = f"storage.googleapis.com/public-download-files/hgnc/archive/archive/{archive}/json/hgnc_complete_set_{version}.json"
    return url


def mane_release(wildcards: Wildcards) -> str | None:
    """The MANE release whose summary sets the MANE tags of the source, if any."""
    key = get_source_config_key(wildcards)
    return config["sources"][key].get("mane", {}).get("release")


def get_mane_summary(wildcards: Wildcards) -> dict[str, str]:
    release = mane_release(wildcards)
    return {
        "summary": f"results/mane/{release}/MANE.GRCh38.v{release}.summary.txt.gz",
        "versions": f"results/mane/{release}/README_versions.txt",
    }


def cdot_input_mapping(wildcards: Wildcards) -> dict[str, str]:
    alias = get_alias(wildcards)
    if len(transcripts_to_fix_with_nuccore(wildcards)) > 0:
        cds = ".cds"
    else:
        cds = ""
    mane = ".mane" if mane_release(wildcards) else ""
    cdot_files = {
        "cdot": f"results/{alias}/cdot/{alias}.cdot.hgnc{cds}{mane}.json.gz",
    }
    if source_kind(wildcards) == "refseq":
        cdot_files.update(
            **{
                "cdot_mt": f"results/{alias}/cdot/{alias}-from-ensembl.chrMT.json",
                "cdot_grafted": f"results/{alias}/cdot/{alias}.cdot.grafted.json.gz",
            }
        )
    return cdot_files


def ensembl_cdot(wildcards: Wildcards) -> str:
    alias_ensembl = used_alias(wildcards, "ensembl")
    return f"results/{alias_ensembl}/cdot/{alias_ensembl}.cdot.hgnc.json.gz"


def get_mehari_input(wildcards: Wildcards) -> dict[str, str]:
    alias = get_alias(wildcards)
    result = {
        "seqrepo_instance": f"results/{alias}/seqrepo/master",
        **cdot_input_mapping(wildcards),
    }
    if wildcards.assembly == "GRCh37":
        other = newest_alias("GRCh38", source_kind(wildcards))
        result.update({"tags": f"results/{other}/cdot/{other}.tags.tsv"})
    return result


def get_mehari_cdot_param_string(wildcards: Wildcards, input: InputFiles) -> str:
    cdot_files = cdot_input_mapping(wildcards)
    expected_input_keys = input.keys()
    for key in cdot_files.keys():
        assert key in expected_input_keys
    return " ".join(f"--annotation {path}" for path in cdot_files.values())


def get_mehari_check_cdot_param_string(wildcards: Wildcards, input: InputFiles) -> str:
    cdot_files = cdot_input_mapping(wildcards)
    expected_input_keys = input.keys()
    for key in cdot_files.keys():
        assert key in expected_input_keys
    return " ".join(f"--cdot {path}" for path in cdot_files.values())


def transcripts_to_fix_with_nuccore(wildcards: Wildcards) -> list[str]:
    key = get_source_config_key(wildcards)
    return list(sorted(set(config["sources"][key].get("fixes", {}).get("fix_cds", []))))


def transcripts_to_lookup_ensembl_ids_for(wildcards: Wildcards) -> list[str]:
    key = get_source_config_key(wildcards)
    return list(
        sorted(set(config["sources"][key].get("fixes", {}).get("add_from_ensembl", [])))
    )


def transcripts_to_fix_polyA(wildcards: Wildcards) -> list[str]:
    key = get_source_config_key(wildcards)
    return list(sorted(set(config["sources"][key].get("fixes", {}).get("polyA", []))))


def parse_missing_accessions(_wildcards, input) -> list[str]:
    with open(input.missing_accessions) as file:
        accessions = {s.strip() for s in file}
    return list(sorted(accessions))


def get_genes_to_disease_download_url() -> str:
    url = "https://github.com/obophenotype/human-phenotype-ontology/releases/download/{release}/genes_to_disease.txt"
    release = config["human-phenotype-ontology"]["genes_to_disease"]["release"]
    return url.format(release=release)


def get_known_issues(wildcards: Wildcards) -> list[str]:
    key = get_source_config_key(wildcards)
    return config["sources"][key].get("known_issues", [])


def get_mehari_docker_url() -> str:
    return config["mehari"]["docker"]
