rule mehari_build_txs_db:
    input:
        unpack(get_mehari_input),
    output:
        txs="results/mehari/{alias}/{seqrepo}/txs.bin.zst",
        report="results/mehari/{alias}/{seqrepo}/txs.bin.zst.report.jsonl",
    params:
        mane=lambda wildcards, input: (
            f"--path-mane-txs-tsv {input.mane_txs}"
            if wildcards.alias == "GRCh37"
            else ""
        ),
        cdot=lambda wildcards, input: (
            f"--path-cdot-json {input.cdot} --path-cdot-json {input.cdot_hgnc}"
            if config["hgnc"]["cdot-mode"] == "create"
            else f"--path-cdot-json {input.cdot_hgnc}"
        ),
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
        --path-cdot-json {input.cdot_mt} \
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
