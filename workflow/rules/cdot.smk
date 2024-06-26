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
        "../scripts/cdot/extract_tags.py"


rule fetch_incorrect_entries:
    params:
        transcripts=transcripts_to_fix_start_stop_codons_for,
    output:
        xml="results/for-fix/{alias}/nuccore.xml.gz",
    log:
        "logs/{alias}/transcripts/fetch_incorrect_entries.log",
    script:
        "../scripts/cdot/fetch_incorrect_entries.py"


rule fix_incorrect_entries:
    input:
        xml="results/for-fix/{alias}/nuccore.xml.gz",
        cdot="results/transcripts/cdot/{alias}.json.gz",
    output:
        cdot="results/for-fix/{alias}/cdot.json.gz",
        report="results/report/{alias}/fix_incorrect_entries.tsv",
    log:
        "logs/{alias}/transcripts/fix_incorrect_entries.log",
    script:
        "../scripts/cdot/fix_incorrect_entries.py"


rule cdot_from_hgnc_complete_set:
    input:
        cdot="results/for-fix/{alias}/cdot.json.gz",
        hgnc="results/hgnc/hgnc_complete_set.json",
    output:
        cdot="results/transcripts/cdot/{alias}.hgnc.json.gz",
        report=report(
            "results/report/{alias}/cdot_hgnc_update.tsv",
            category="{alias}",
        ),
    log:
        "logs/{alias}/transcripts/update_cdot_with_hgnc_complete_set.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/update_from_hgnc.py"


rule cdot_chrMT:
    input:
        cdot="results/transcripts/cdot/{alias}-ensembl.json.gz",
    output:
        cdot="results/transcripts/cdot/{alias}-ensembl.chrMT.json",
    log:
        "logs/{alias}-ensembl/transcripts/extract_chrMT.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/extract_chrmt.py"


rule lookup_ensembl_ids_for_refseq_ids:
    params:
        refseq_ids=transcripts_to_lookup_ensembl_ids_for,
    output:
        tsv="results/for-fix/{alias}/refseq_id_to_ensembl_id.tsv",
    log:
        "logs/{alias}-ensembl/transcripts/lookup_ensembl_ids_for_refseq_ids.log",
    script:
        "../scripts/lookup_ensembl_ids_for_refseq_ids.py"


rule cdot_graft_ensembl_ids_for_certain_refseq_ids:
    input:
        cdot="results/transcripts/cdot/{alias}-ensembl.json.gz",
        lookup="results/for-fix/{alias}/refseq_id_to_ensembl_id.tsv",
    output:
        cdot="results/transcripts/cdot/{alias}-ensembl.grafted.json.gz",
    log:
        "logs/{alias}-ensembl/transcripts/graft_ensembl_ids.log",
    script:
        "../scripts/cdot/graft_ensembl_ids.py"
