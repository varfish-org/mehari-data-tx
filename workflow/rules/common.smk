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
    if alias == "GRCh38":
        return {
            # CASP8AP2
            "NM_001137667.2",
            "NM_001137668.2",
            "NM_012115.4",
            # DAG1
            "NM_001177639.3",
            # FOXO6
            "NM_001291281.3",
            # POLK
            "NM_001345921.3",
            "NM_001345922.3",
            "NM_001372044.2",
            "NM_001387110.3",
            "NM_001387111.3",
            "NM_001395899.1",
            "NM_001395901.1",
            "NM_001395902.1",
            "NM_016218.6",
            # SAMD1
            "NM_138352.3",
            # MUC19
            "NM_173600.2",
            # older versions in cdot than available in the current cdna ref releases:
            "NM_001401501.1",
            "NM_002457.4",
            "NM_003890.2",
            "NM_012234.6",
            "NM_017940.6",
            # Invalid CDS length in cdot:
            # Most (or all?) of these are marked as CDS join(a..b,c..d) in nuccore
            "NM_001172437.2",
            "NM_001184961.1",
            "NM_015068.3",
            "NM_182705.2",
            "NM_001145051.2",
            "NM_001301020.1",
            "NM_004152.3",
            "NM_001301302.1",
            "NM_002537.3",
            "NM_001134939.1",
            "NM_001301371.1",
            "NM_016178.2",
            "NM_053005.5",
        }
    else:
        # TODO: transcript retrieval for GRCh37, uses a dummy value atm
        return {"NM_182705.2"}
