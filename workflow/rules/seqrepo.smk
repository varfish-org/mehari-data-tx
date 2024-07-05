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
        # or `.value.reason | contains("MissingSequence")`, but if there are more reasons a transcript has been discarded,
        # it makes no sense to fetch transcripts with missing sequence if they get discarded for a different reason anyway
        """(jq -r 'select(.value.reason == "MissingSequence").value.id.value' {input.txs_db_report} > {output.missing_txt}) >{log} 2>&1"""


rule fetch_missing_sequence:
    output:
        missing_fasta=ensure(
            "results/{assembly}-{source}/mehari/seqrepo/missing/{accession}.fasta",
            non_empty=True,
        ),
    resources:
        ratelimit=1,
    params:
        accession=lambda wildcards: wildcards.accession,
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
        checksum="results/{assembly}-{source}/mehari/seqrepo/missing.fasta.md5",
        new_checksum=temp(
            "results/{assembly}-{source}/mehari/seqrepo/missing.fasta.md5.new.tmp"
        ),
    log:
        "logs/{assembly}-{source}/mehari/seqrepo/fetch-missing.log",
    run:
        if os.path.exists(output.checksum):
            old_checksum = open(output.checksum).readline().strip().split()[0]
        else:
            old_checksum = None
        if len(input) > 0:
            shell("cat {input} > {output.missing_fasta}")
            shell("md5sum {output.missing_fasta} > {output.new_checksum}")
            new_checksum = open(output.new_checksum).readline().strip().split()[0]
        else:
            shell("touch {output.missing_fasta}")
            shell("md5sum {output.missing_fasta} > {output.checksum}")
            shell("touch {output.new_checksum}")
            new_checksum = None
        if old_checksum != new_checksum:
            shell("cp -up {output.new_checksum} {output.checksum}")


rule fix_missing_sequences_in_seqrepo:
    input:
        seqrepo_root="results/{assembly}-{source}/seqrepo",
        seqrepo_instance="results/{assembly}-{source}/seqrepo/master",
        missing_fasta_checksum="results/{assembly}-{source}/mehari/seqrepo/missing.fasta.md5",
    output:
        seqrepo_root_fixed=directory("results/{assembly}-{source}/seqrepo_fixed"),
        seqrepo_fixed_instance=directory(
            "results/{assembly}-{source}/seqrepo_fixed/master"
        ),
    params:
        missing_fasta=lambda wildcards, input: input.missing_fasta_checksum.rstrip(
            ".md5"
        ),
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
          load --instance-name master --namespace {params.refseq_namespace} {params.missing_fasta} | tail
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root_fixed} \
          load --instance-name master --namespace {params.ensembl_namespace} {params.missing_fasta} | tail
        """
