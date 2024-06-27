rule get_cdot_transcripts:
    output:
        "results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.json.gz",
    params:
        url=get_cdot_download_url,
    conda:
        "../envs/base.yaml"
    log:
        "logs/{assembly}-{source}/get_cdot_transcripts.log",
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
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.json.gz",
    output:
        mane_txs="results/{assembly}-{source}/cdot/{assembly}-{source}.mane-txs.tsv",
    log:
        "logs/{assembly}-{source}/cdot/mane_txs_for_grch37.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/extract_tags.py"


rule fetch_incorrect_entries:
    params:
        transcripts=transcripts_to_fix_start_stop_codons_for,
    output:
        xml="results/{assembly}-{source}/cdot/nuccore.xml.gz",
    log:
        "logs/{assembly}-{source}/cdot/fetch_incorrect_entries.log",
    script:
        "../scripts/cdot/fetch_incorrect_entries.py"


rule fix_incorrect_entries:
    input:
        xml="results/{assembly}-{source}/cdot/nuccore.xml.gz",
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.json.gz",
    output:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.fixed.json.gz",
        report="results/{assembly}-{source}/report/fix_incorrect_entries.tsv",
    log:
        "logs/{assembly}-{source}/cdot/fix_incorrect_entries.log",
    script:
        "../scripts/cdot/fix_incorrect_entries.py"


rule cdot_from_hgnc_complete_set:
    input:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.fixed.json.gz",
        hgnc="results/hgnc/hgnc_complete_set.json",
    output:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.fixed.hgnc.json.gz",
        report=report(
            "results/report/{assembly}-{source}/cdot_hgnc_update.tsv",
            category="{assembly}-{source}",
        ),
    log:
        "logs/{assembly}-{source}/cdot/update_cdot_with_hgnc_complete_set.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/update_from_hgnc.py"


rule cdot_chrMT:
    input:
        cdot=ensembl_cdot,
    output:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}-from-ensembl.chrMT.json",
    log:
        "logs/{assembly}-{source}/cdot/extract_chrMT.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/extract_chrmt.py"


rule lookup_ensembl_ids_for_refseq_ids:
    params:
        refseq_ids=transcripts_to_lookup_ensembl_ids_for,
    output:
        tsv="results/{assembly}-{source}/lookup/refseq_id_to_ensembl_id.tsv",
    log:
        "logs/{assembly}-{source}/lookup/lookup_ensembl_ids_for_refseq_ids.log",
    script:
        "../scripts/lookup_ensembl_ids_for_refseq_ids.py"


rule cdot_graft_ensembl_ids_for_certain_refseq_ids:
    input:
        cdot=ensembl_cdot,
        lookup="results/{assembly}-{source}/lookup/refseq_id_to_ensembl_id.tsv",
    output:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}-from-ensembl.grafted.json.gz",
    log:
        "logs/{assembly}-{source}/cdot/graft_ensembl_ids.log",
    script:
        "../scripts/cdot/graft_ensembl_ids.py"
