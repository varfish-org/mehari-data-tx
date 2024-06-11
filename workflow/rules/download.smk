rule get_ensembl_sequence:
    output:
        "results/references/ensembl/{alias}.fasta",
    params:
        species=get_ensembl_sequence_param("species"),
        build=get_ensembl_sequence_param("build"),
        release=get_ensembl_sequence_param("release"),
        datatype=get_ensembl_sequence_param("datatype"),
    log:
        "logs/{alias}/get_ensembl_sequence.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v3.11.0/bio/reference/ensembl-sequence"


rule get_refseq_sequence:
    output:
        refseq_seq="results/references/refseq/{alias}.fasta.gz",
    params:
        species=lambda wildcards: config["reference"][wildcards.alias]["refseq"][
            "species"
        ],
        species_name=lambda wildcards: config["reference"][wildcards.alias]["refseq"][
            "species_name"
        ],
    log:
        "logs/{alias}/get_refseq_sequence.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/download_refseq.py"


rule get_cdot_transcripts:
    output:
        "results/transcripts/cdot/{alias}.json.gz",
    params:
        url=get_cdot_download_url,
    conda:
        "../envs/base.yaml"
    log:
        "logs/{alias}/get_cdot_transcripts.log",
    shell:
        """wget --quiet -O {output} {params.url} 2> {log}"""


rule get_hgnc_complete_set:
    output:
        "results/hgnc/hgnc_complete_set.json",
    params:
        url=get_hgnc_complete_set_download_url,
    conda:
        "../envs/base.yaml"
    log:
        "logs/hgnc/get_hgnc_complete_set.log",
    shell:
        """wget --quiet -O {output} {params.url} 2> {log}"""
