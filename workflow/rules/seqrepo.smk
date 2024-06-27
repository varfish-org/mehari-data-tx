rule initialize_seqrepo:
    input:
        ensembl="results/{assembly}-ensembl/reference/{assembly}-ensembl.fasta",
        refseq="results/{assembly}-refseq/reference/{assembly}-refseq.fasta.gz",
    output:
        seqrepo_root=directory("results/{assembly}-{source}/seqrepo"),
        seqrepo_instance=directory("results/{assembly}-{source}/seqrepo/master"),
    params:
        ensembl_namespace=config["namespaces"]["ensembl"],
        refseq_namespace=config["namespaces"]["refseq"],
    log:
        "logs/{assembly}-{source}/seqrepo/initialize.log",
    conda:
        "../envs/seqrepo.yaml"
    shell:
        """
        seqrepo --root-directory {output.seqrepo_root} init --instance-name master

        # seqrepo load is too verbose and we cannot silence it
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root} \
            load --instance-name master --namespace {params.ensembl_namespace} \
            {input.ensembl} \
        | tail

        # seqrepo load is too verbose and we cannot silence it
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root} \
            load --instance-name master --namespace {params.refseq_namespace} \
            {input.refseq} \
        | tail
        """


checkpoint detect_missing_sequences:
    input:
        txs_db_report="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.report.jsonl",
    output:
        missing_txt=report(
            "results/{assembly}-{source}/mehari/seqrepo/missing.txt",
            category="{assembly}-{source}",
            subcategory="seqrepo",
        ),
    log:
        "logs/{assembly}-{source}/seqrepo/detect-missing.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """(jq -r 'select(.value.reason == "MissingSequence").value.id.value' {input.txs_db_report} > {output.missing_txt}) >{log} 2>&1"""


rule fetch_missing_sequence:
    output:
        missing_fasta="results/{assembly}-{source}/mehari/seqrepo/missing/{accession}.fasta",
    params:
        accession=lambda wildcards: wildcards.accession,
    resources:
        ratelimit=1,
    cache: "omit-software"
    log:
        "logs/{assembly}-{source}/mehari/seqrepo/fetch-missing/{accession}.log",
    conda:
        "../envs/seqrepo.yaml"
    script:
        "../scripts/fetch_missing_sequence.py"


rule aggregate_missing_sequences:
    input:
        missing_sequence_files,
    output:
        missing_fasta="results/{assembly}-{source}/mehari/seqrepo/missing.fasta",
    log:
        "logs/{assembly}-{source}/mehari/seqrepo/fetch-missing.log",
    shell:
        """
        (
        cat {input} > {output.missing_fasta}
        ) >{log} 2>&1
        """


rule fix_missing_sequences_in_seqrepo:
    input:
        seqrepo_root="results/{assembly}-{source}/seqrepo",
        seqrepo_instance="results/{assembly}-{source}/seqrepo/master",
        missing_fasta="results/{assembly}-{source}/mehari/seqrepo/missing.fasta",
    output:
        seqrepo_root_fixed=directory("results/{assembly}-{source}/seqrepo_fixed"),
        seqrepo_fixed_instance=directory(
            "results/{assembly}-{source}/seqrepo_fixed/master"
        ),
    params:
        refseq_namespace=config["namespaces"]["refseq"],
        ensembl_namespace=config["namespaces"]["ensembl"],
    conda:
        "../envs/seqrepo.yaml"
    log:
        "logs/{assembly}-{source}/seqrepo/fix-missing.log",
    shell:
        """
        cp -r {input.seqrepo_root}/* {output.seqrepo_root_fixed}
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root_fixed} \
          load --instance-name master --namespace {params.refseq_namespace} {input.missing_fasta} | tail
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root_fixed} \
          load --instance-name master --namespace {params.ensembl_namespace} {input.missing_fasta} | tail
        """
