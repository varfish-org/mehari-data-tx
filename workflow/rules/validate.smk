rule known_issues_tsv:
    output:
        known_issues="results/{assembly}-{source}/fixes/known_issues.tsv",
    params:
        known_issues=get_known_issues,
    conda:
        "../envs/datastuff.yaml"
    log:
        "logs/{assembly}-{source}/fixes/known_issues.log",
    script:
        "../scripts/known_issues.py"


rule check_mehari_db:
    input:
        unpack(cdot_input_mapping),
        tx_db="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst",
        tx_db_report="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.report.jsonl",
        tx_db_sha256="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.sha256",
        tx_db_report_sha256="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.report.jsonl.sha256",
        hgnc="results/hgnc/hgnc_complete_set.json",
        genes_to_disease="results/human-phenotype-ontology/genes_to_disease_with_hgnc_id.tsv",
        known_issues="results/{assembly}-{source}/fixes/known_issues.tsv",
    output:
        # stats="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.stats.tsv",
        report=report(
            "results/{assembly}-{source}/mehari/seqrepo/report/mehari_db_check.txt",
            category="{assembly}-{source}",
        ),
    params:
        cdot=get_mehari_check_cdot_param_string,
    log:
        "logs/{assembly}-{source}/mehari/seqrepo/check.log",
    container:
        get_mehari_docker_url()
    shell:
        """(
        mehari db check \
        --db {input.tx_db} \
        --hgnc {input.hgnc} \
        --disease-genes {input.genes_to_disease} \
        {params.cdot} \
        --known-issues {input.known_issues} \
        --output {output.report} \
        ) >{log} 2>&1"""


rule generate_table_with_discards_to_review:
    input:
        discarded="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.report.jsonl",
    output:
        table="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.discards.tsv",
    log:
        "logs/{assembly}-{source}/mehari/seqrepo/discards_table.log",
    conda:
        "../envs/datastuff.yaml"
    shell:
        """(
        echo -e "idtype\tid\tgene_name\treason\ttags" > {output.table}
        jq -r 'select( .type == "Discard" ) | [.value.id.type, .value.id.value, .value.gene_name, .value.reason, (.value.tags // [] | map(.Other? // .) | join(","))] | @tsv' {input.discarded} >> {output.table}
        ) >{log} 2>&1"""


rule datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/report.datavzrd.yaml"),
        mehari_check_db_stats="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.stats.tsv",
        fix_incorrect_cds="results/{assembly}-{source}/cdot/fix_incorrect_cds.tsv",
        refseq_id_to_ensembl_id="results/{assembly}-{source}/lookup/refseq_id_to_ensembl_id.tsv",
        discards_of_interest="results/{assembly}-{source}/mehari/seqrepo/txs.bin.zst.discards.tsv",
    params:
        extra="",
    output:
        report(
            directory("results/{assembly}-{source}/datavzrd-report/seqrepo"),
            htmlindex="index.html",
            category="{assembly}-{source}",
        ),
    log:
        "logs/datavzrd_report/{assembly}-{source}/seqrepo/check_mehari_db.log",
    wrapper:
        "v5.0.2/utils/datavzrd"
