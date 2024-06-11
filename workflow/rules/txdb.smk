rule cdot_chrMT:
    input:
        cdot="results/transcripts/cdot/GRCh38-ensembl.json.gz",
    output:
        cdot="results/transcripts/cdot/GRCh38-ensembl.chrMT.json",
    log:
        "logs/GRCh38-ensembl/transcripts/extract_chrMT.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot_extract_chrmt.py"


rule mane_txs_for_grch37:
    input:
        cdot="results/transcripts/cdot/GRCh37.json.gz",
    output:
        mane_txs="results/transcripts/cdot/GRCh37/mane-txs.tsv",
    log:
        "logs/GRCh37/transcripts/mane_txs_for_grch37.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot_json_to_tags.py"


rule update_cdot_with_hgnc_complete_set:
    input:
        cdot="results/transcripts/cdot/{alias}.json.gz",
        hgnc="results/hgnc/hgnc_complete_set.json",
    output:
        cdot="results/transcripts/cdot/{alias}.hgnc.json.gz",
    log:
        "logs/{alias}/transcripts/update_cdot_with_hgnc_complete_set.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot_from_hgnc.py"


rule mehari_build_txs_db:
    input:
        unpack(get_mehari_input),
    output:
        txs="results/mehari/{alias}/{seqrepo}/txs.bin.zst",
        report="results/mehari/{alias}/{seqrepo}/txs.bin.zst.report.jsonl",
    params:
        extra=lambda wildcards, input: (
            f"--path-mane-txs-tsv {input.mane_txs}"
            if wildcards.alias == "GRCh37"
            else ""
        ),
        genome_release=lambda wildcards: wildcards.alias.lower(),
    log:
        "logs/{alias}/mehari/{seqrepo}/build_txs_db.log",
    conda:
        "../envs/mehari.yaml"
    shell:
        """
        mehari db create \
        --path-out {output.txs} \
        --path-seqrepo-instance {input.seqrepo_instance} \
        --path-cdot-json {input.cdot} \
        --path-cdot-json {input.cdot_hgnc} \
        --path-cdot-json {input.cdot_mt} \
        --genome-release {params.genome_release} \
        {params.extra} 2> {log}
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
