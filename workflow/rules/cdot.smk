rule get_cdot_transcripts:
    output:
        "results/{alias}/cdot/{alias}.cdot.json.gz",
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
        cdot="results/{alias}/cdot/{alias}.cdot.json.gz",
    output:
        mane_txs="results/{alias}/cdot/{alias}.mane-txs.tsv",
    log:
        "logs/{alias}/cdot/mane_txs_for_grch37.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/extract_tags.py"


rule fetch_incorrect_entries:
    params:
        transcripts=transcripts_to_fix_start_stop_codons_for,
    output:
        xml="results/{alias}/cdot/nuccore.xml.gz",
    log:
        "logs/{alias}/cdot/fetch_incorrect_entries.log",
    script:
        "../scripts/cdot/fetch_incorrect_entries.py"


rule fix_incorrect_entries:
    input:
        xml="results/{alias}/cdot/nuccore.xml.gz",
        cdot="results/{alias}/cdot/{alias}.cdot.json.gz",
    output:
        cdot="results/{alias}/cdot/{alias}.cdot.fixed.json.gz",
        report="results/{alias}/report/fix_incorrect_entries.tsv",
    log:
        "logs/{alias}/cdot/fix_incorrect_entries.log",
    script:
        "../scripts/cdot/fix_incorrect_entries.py"


rule cdot_from_hgnc_complete_set:
    input:
        cdot="results/{alias}/cdot/{alias}.cdot.fixed.json.gz",
        hgnc="results/hgnc/hgnc_complete_set.json",
    output:
        cdot="results/{alias}/cdot/{alias}.cdot.fixed.hgnc.json.gz",
        report=report(
            "results/report/{alias}/cdot_hgnc_update.tsv",
            category="{alias}",
        ),
    log:
        "logs/{alias}/cdot/update_cdot_with_hgnc_complete_set.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/update_from_hgnc.py"


rule cdot_chrMT:
    input:
        cdot=ensembl_cdot,
    output:
        cdot="results/{alias}/cdot/{alias}-from-ensembl.chrMT.json",
    log:
        "logs/{alias}/cdot/extract_chrMT.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/extract_chrmt.py"


rule lookup_ensembl_ids_for_refseq_ids:
    params:
        refseq_ids=transcripts_to_lookup_ensembl_ids_for,
    output:
        tsv="results/{alias}/lookup/refseq_id_to_ensembl_id.tsv",
    log:
        "logs/{alias}/lookup/lookup_ensembl_ids_for_refseq_ids.log",
    script:
        "../scripts/lookup_ensembl_ids_for_refseq_ids.py"


rule cdot_graft_ensembl_ids_for_certain_refseq_ids:
    input:
        cdot=ensembl_cdot,
        lookup="results/{alias}/lookup/refseq_id_to_ensembl_id.tsv",
    output:
        cdot="results/{alias}/cdot/{alias}-from-ensembl.grafted.json.gz",
    log:
        "logs/{alias}/cdot/graft_ensembl_ids.log",
    script:
        "../scripts/cdot/graft_ensembl_ids.py"
