rule mehari_build_txs_db:
    input:
        unpack(get_mehari_input),
    output:
        txs="results/{alias}/mehari/{seqrepo}/txs.bin.zst",
        report="results/{alias}/mehari/{seqrepo}/txs.bin.zst.report.jsonl",
    params:
        mane=lambda wildcards, input: (
            f"--path-mane-txs-tsv {input.mane_txs}"
            if wildcards.alias == "GRCh37"
            else ""
        ),
        cdot=get_mehari_cdot_param_string,
        genome_release=lambda wildcards: wildcards.alias.lower(),
    log:
        "logs/{alias}/mehari/{seqrepo}/build_txs_db.log",
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
