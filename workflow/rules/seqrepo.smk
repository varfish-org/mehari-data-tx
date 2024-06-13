rule initialize_seqrepo:
    input:
        ensembl="results/references/ensembl/{alias}.fasta",
        refseq="results/references/refseq/{alias}.fasta.gz",
    output:
        seqrepo_root=directory("results/seqrepo/{alias}"),
        seqrepo_instance=directory("results/seqrepo/{alias}/master"),
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


rule fix_missing:
    input:
        seqrepo_root="results/seqrepo/{alias}",
        seqrepo_instance="results/seqrepo/{alias}/master",
        txs_db="results/mehari/{alias}/seqrepo/txs.bin.zst",
        txs_db_report="results/mehari/{alias}/seqrepo/txs.bin.zst.report.jsonl",
    output:
        seqrepo_root_fixed=directory("results/seqrepo_fixed/{alias}"),
        seqrepo_fixed_instance=directory("results/seqrepo_fixed/{alias}/master"),
        missing_txt="results/for-fix/{alias}/missing.txt",
        missing_fasta="results/for-fix/{alias}/missing.fasta",
    params:
        refseq_namespace=config["namespaces"]["refseq"],
    conda:
        "../envs/seqrepo.yaml"
    log:
        "logs/{alias}/seqrepo/fix-missing.log",
    shell:
        """
        set -euo pipefail
        IFS=$'\n\t'
        set -x
        jq -r 'select(.reason == "MissingSequence").id' {input.txs_db_report} > {output.missing_txt}
        cp -r {input.seqrepo_root}/* {output.seqrepo_root_fixed}
        2>{log} xargs -a {output.missing_txt} \
         seqrepo --root-directory {output.seqrepo_root_fixed} \
          fetch-load --instance-name master --namespace {params.refseq_namespace} >> {output.missing_fasta}
        """
