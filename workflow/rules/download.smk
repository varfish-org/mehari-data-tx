rule get_ensembl_sequence:
    output:
        "results/{assembly}-{source}/reference/{assembly}-{source}.fasta",
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


rule get_refseq_installed_files:
    output:
        "results/{assembly}-{source}/reference/files.installed",
    params:
        species=get_refseq_sequence_param("species"),
        species_name=get_refseq_sequence_param("species_name"),
    log:
        "logs/{assembly}-{source}/get_refseq_installed_files.log",
    conda:
        "../envs/base.yaml"
    shell:
        """(
        curl --silent "https://ftp.ncbi.nih.gov/refseq/{params.species}/mRNA_Prot/{params.species_name}.files.installed" > {output}
        ) >{log} 2>&1"""


rule get_refseq_sequence:
    input:
        refseq_installed="results/{assembly}-{source}/reference/files.installed",
    output:
        refseq_seq="results/{assembly}-{source}/reference/{assembly}-{source}.fasta.gz",
    params:
        species=get_refseq_sequence_param("species"),
        species_name=get_refseq_sequence_param("species_name"),
    cache: "omit-software"
    log:
        "logs/{assembly}-{source}/get_refseq_sequence.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/download_refseq.py"
