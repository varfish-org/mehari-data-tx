rule mehari_build_txs_db:
    input:
        unpack(get_mehari_input),
    output:
        txs="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst",
        report=report(
            "results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.report.jsonl",
            category="{assembly}-{source}",
        ),
    threads: workflow.cores / 2
    params:
        mane=lambda wildcards, input: (
            f"--path-mane-txs-tsv {input.tags}"
            if wildcards.assembly == "GRCh37"
            else ""
        ),
        cdot=get_mehari_cdot_param_string,
        assembly=genome_assembly,
        assembly_version=genome_assembly_version_parameter,
        cdot_version=cdot_version,
        transcript_source=transcript_source,
        transcript_source_version=transcript_source_version_parameter,
    log:
        "logs/{assembly}-{source}/mehari/seqrepo/build_txs_db.log",
    benchmark:
        "benchmarks/{assembly}-{source}/mehari/seqrepo/build_txs_db.tsv"
    container:
        get_mehari_docker_url()
    shell:
        """
        mehari db create \
        --path-out {output.txs} \
        --path-seqrepo-instance {input.seqrepo_instance} \
        {params.cdot} \
        --cdot-version {params.cdot_version} \
        --assembly {params.assembly} \
        {params.assembly_version} \
        --transcript-source {params.transcript_source} \
        {params.transcript_source_version} \
        --threads {threads} \
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
