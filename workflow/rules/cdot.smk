rule get_cdot_transcripts:
    output:
        "results/transcripts/cdot/{alias}.json.gz",
    params:
        url=get_cdot_download_url,
    conda:
        "../envs/base.yaml"
    log:
        "logs/{alias}/get_cdot_transcripts.log",
    shell:
        """wget --quiet -O {output} {params.url} 2> {log}"""


rule get_hgnc_complete_set:
    output:
        "results/hgnc/hgnc_complete_set.json",
    params:
        url=get_hgnc_complete_set_download_url,
    conda:
        "../envs/base.yaml"
    log:
        "logs/hgnc/get_hgnc_complete_set.log",
    shell:
        """wget --quiet -O {output} {params.url} 2> {log}"""


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


rule fetch_incorrect_entries:
    params:
        transcripts=transcripts_to_fix_start_stop_codons_for,
    output:
        xml="results/for-fix/{alias}/nuccore.xml.gz",
    script:
        "../scripts/cdot_fetch_incorrect_entries.py"


rule fix_incorrect_entries:
    input:
        xml="results/for-fix/{alias}/nuccore.xml.gz",
        cdot="results/transcripts/cdot/{alias}.json.gz",
    output:
        cdot="results/for-fix/{alias}/cdot.json.gz",
    script:
        "../scripts/cdot_fix_incorrect_entries.py"


rule cdot_from_hgnc_complete_set:
    input:
        cdot="results/for-fix/{alias}/cdot.json.gz",
        hgnc="results/hgnc/hgnc_complete_set.json",
    output:
        cdot="results/transcripts/cdot/{alias}.hgnc.json.gz",
    params:
        mode=config["hgnc"]["cdot-mode"],
    log:
        "logs/{alias}/transcripts/update_cdot_with_hgnc_complete_set.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot_from_hgnc.py"


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
