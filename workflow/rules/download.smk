rule get_ensembl_sequence:
    output:
        "results/references/ensembl/{assembly}-{source}.fasta",
    params:
        species=get_ensembl_sequence_param("species"),
        build=get_ensembl_sequence_param("build"),
        release=get_ensembl_sequence_param("release"),
        datatype=get_ensembl_sequence_param("datatype"),
    log:
        "logs/{assembly}-{source}/get_ensembl_sequence.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v3.11.0/bio/reference/ensembl-sequence"


rule get_refseq_sequence:
    output:
        refseq_seq="results/references/refseq/{assembly}-{source}.fasta.gz",
    params:
        species=get_refseq_sequence_param("species"),
        species_name=get_refseq_sequence_param("species_name"),
    log:
        "logs/{assembly}-{source}/get_refseq_sequence.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/download_refseq.py"
