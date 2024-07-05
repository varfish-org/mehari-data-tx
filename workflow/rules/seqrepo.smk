import os


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


rule detect_missing_sequences:
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
        # or `.value.reason | contains("MissingSequence")`, but if there are more reasons a transcript has been discarded,
        # it makes no sense to fetch transcripts with missing sequence if they get discarded for a different reason anyway
        """(jq -r 'select(.value.reason == "MissingSequence").value.id.value' {input.txs_db_report} > {output.missing_txt}) >{log} 2>&1"""


rule aggregate_missing_sequences:
    input:
        missing_accessions="results/{assembly}-{source}/mehari/seqrepo/missing.txt",
    output:
        missing_fasta="results/{assembly}-{source}/mehari/seqrepo/missing.fasta",
        missing_fasta_fai="results/{assembly}-{source}/mehari/seqrepo/missing.fasta.fai",
    params:
        missing_accessions=parse_missing_accessions,
    resources:
        ratelimit=1,
    log:
        "logs/{assembly}-{source}/mehari/seqrepo/fetch-missing.log",
    conda:
        "../envs/efetch.yaml"
    script:
        "../scripts/fetch_missing_sequences.py"


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
