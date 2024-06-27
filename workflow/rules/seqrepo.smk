rule initialize_seqrepo:
    input:
        ensembl="results/references/ensembl/{alias}.fasta",
        refseq="results/references/refseq/{alias}.fasta.gz",
    output:
        seqrepo_root=directory("results/{alias}/seqrepo"),
        seqrepo_instance=directory("results/{alias}/seqrepo/master"),
    params:
        ensembl_namespace=config["namespaces"]["ensembl"],
        refseq_namespace=config["namespaces"]["refseq"],
    log:
        "logs/{alias}/seqrepo/initialize.log",
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
        txs_db_report="results/{alias}/mehari/seqrepo/txs.bin.zst.report.jsonl",
    output:
        missing_txt=report(
            "results/{alias}/mehari/seqrepo/missing.txt",
            category="{alias}",
            subcategory="seqrepo",
        ),
    log:
        "logs/{alias}/seqrepo/detect-missing.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """(jq -r 'select(.value.reason == "MissingSequence").value.id.value' {input.txs_db_report} > {output.missing_txt}) >{log} 2>&1"""


rule fetch_missing_sequence:
    input:
        missing_txt="results/{alias}/mehari/seqrepo/missing/{accession}.txt",
    output:
        missing_fasta="results/{alias}/mehari/seqrepo/missing/{accession}.fasta",
    resources:
        ratelimit=1,
    cache: "omit-software"
    log:
        "logs/{alias}/mehari/seqrepo/fetch-missing/{accession}.log",
    conda:
        "../envs/seqrepo.yaml"
    script:
        "../scripts/fetch_missing_sequences.py"


rule fetch_missing_sequences:
    input:
        missing_sequence_files,
    output:
        missing_fasta="results/{alias}/mehari/seqrepo/missing.fasta",
    log:
        "logs/{alias}/mehari/seqrepo/fetch-missing.log",
    shell:
        """
        (
        cat {input} > {output.missing_fasta}
        ) >{log} 2>&1
        """


rule fix_missing_sequences_in_seqrepo:
    input:
        seqrepo_root="results/{alias}/seqrepo",
        seqrepo_instance="results/{alias}/seqrepo/master",
        missing_fasta="results/{alias}/mehari/seqrepo/missing.fasta",
    output:
        seqrepo_root_fixed=directory("results/{alias}/seqrepo_fixed"),
        seqrepo_fixed_instance=directory("results/{alias}/seqrepo_fixed/master"),
    params:
        refseq_namespace=config["namespaces"]["refseq"],
        ensembl_namespace=config["namespaces"]["ensembl"],
    conda:
        "../envs/seqrepo.yaml"
    log:
        "logs/{alias}/seqrepo/fix-missing.log",
    shell:
        """
        cp -r {input.seqrepo_root}/* {output.seqrepo_root_fixed}
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root_fixed} \
          load --instance-name master --namespace {params.refseq_namespace} {input.missing_fasta} | tail
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root_fixed} \
          load --instance-name master --namespace {params.ensembl_namespace} {input.missing_fasta} | tail
        """
