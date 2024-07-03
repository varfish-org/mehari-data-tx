rule mehari_build_txs_db:
    input:
        unpack(get_mehari_input),
    output:
        txs="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst",
        report=report(
            "results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.report.jsonl",
            category="{assembly}-{source}",
            subcategory="{seqrepo}",
        ),
    threads: workflow.cores
    params:
        mane=lambda wildcards, input: (
            f"--path-mane-txs-tsv {input.tags}"
            if wildcards.assembly == "GRCh37"
            else ""
        ),
        cdot=get_mehari_cdot_param_string,
        genome_release=genome_release,
    log:
        "logs/{assembly}-{source}/mehari/{seqrepo}/build_txs_db.log",
    # conda:
    #     "../envs/mehari.yaml"
    shell:
        """
        mehari db create \
        --path-out {output.txs} \
        --path-seqrepo-instance {input.seqrepo_instance} \
        {params.cdot} \
        --genome-release {params.genome_release} \
        {params.mane} 2> {log}
        """


rule calculate_sha256_sum:
    input:
        "{file}",
    output:
        "{file}.sha256",
    log:
        "logs/sha256sum/{file}.log",
    conda:
        "../envs/checksum.yaml"
    shell:
        "sha256sum {input} > {output} 2> {log}"
