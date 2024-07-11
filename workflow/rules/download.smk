rule get_ensembl_sequence:
    output:
        "results/{assembly}-{source}/reference/{assembly}-{source}.{datatype}.fasta",
    params:
        species=get_ensembl_sequence_param("species"),
        build=get_ensembl_sequence_param("build"),
        release=get_ensembl_sequence_param("release"),
        datatype=lambda wildcards: wildcards.datatype,
    log:
        "logs/{assembly}-{source}/get_ensembl_sequence.{datatype}.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v3.11.0/bio/reference/ensembl-sequence"


rule merge_ensembl_sequence:
    input:
        "results/{assembly}-{source}/reference/{assembly}-{source}.cdna.fasta",
        "results/{assembly}-{source}/reference/{assembly}-{source}.ncrna.fasta",
    output:
        "results/{assembly}-{source}/reference/{assembly}-{source}.fasta",
    log:
        "logs/{assembly}-{source}/merge_sequence.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    conda:
        "../envs/base.yaml"
    shell:
        """
        (
        cat {input} > {output}
        ) >{log} 2>&1
        """


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


rule get_cdot_data:
    output:
        "results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.json.gz",
    params:
        url=get_cdot_download_url,
    conda:
        "../envs/base.yaml"
    log:
        "logs/{assembly}-{source}/get_cdot_transcripts.log",
    shell:
        """curl --silent {params.url} | gzip -dc | jq '.' | gzip -c > {output} 2> {log}"""


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
        """curl --silent {params.url} | jq '.' > {output} 2> {log}"""


rule get_genes_to_disease:
    output:
        "results/human-phenotype-ontology/genes_to_disease.tsv",
    params:
        url=get_genes_to_disease_download_url(),
    conda:
        "../envs/base.yaml"
    log:
        "logs/human-phenotype-ontology/get_genes_to_disease.log",
    shell:
        """wget --quiet -O {output} {params.url} 2> {log}"""


rule add_hgnc_id_to_genes_to_disease:
    input:
        genes_to_disease="results/human-phenotype-ontology/genes_to_disease.tsv",
        hgnc_complete_set="results/hgnc/hgnc_complete_set.json",
    output:
        genes_to_disease="results/human-phenotype-ontology/genes_to_disease_with_hgnc_id.tsv",
    conda:
        "../envs/datastuff.yaml"
    log:
        "logs/human-phenotype-ontology/add_hgnc_id_to_genes_to_disease.log",
    script:
        "../scripts/human_phenotype_ontology/add_hgnc_id_to_genes_to_disease.py"
