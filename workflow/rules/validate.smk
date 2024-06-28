rule dump_mehari_db:
    input:
        tx_db="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst",
    output:
        db_yaml="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.yaml.gz",
    # conda:
    #     "../envs/mehari.yaml"
    log:
        "logs/{assembly}-{source}/mehari/{seqrepo}/dump.log",
    shell:
        """
        (
            mehari db dump --path-db {input.tx_db} | gzip -c > {output.db_yaml}
        ) >{log} 2>&1
        """


rule mehari_transcript_ids:
    input:
        db_yaml="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.yaml.gz",
    output:
        tx_ids="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.transcript_ids.txt",
    log:
        "logs/{assembly}-{source}/mehari/{seqrepo}/transcript_ids.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """yq -r '.txDb.transcripts[].id' <(pigz -dc {input.db_yaml}) > {output.tx_ids} 2> {log}"""


rule cdot_tx_ids:
    input:
        unpack(cdot_input_mapping),
    output:
        tx_ids="results/{assembly}-{source}/cdot/{assembly}-{source}.txids.txt",
    log:
        "logs/{assembly}-{source}/cdot/txids.log",
    shell:
        """
        (
        touch {output.tx_ids}
        for f in {input}; do
            pigz -dcf ${{f}} | jq -r '.transcripts | keys[]' >> {output.tx_ids}
        done
        ) >{log} 2>&1
        """


rule check_mehari_db:
    input:
        kept="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.transcript_ids.txt",
        db_discarded="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.report.jsonl",
        cdot_tx_ids="results/{assembly}-{source}/cdot/{assembly}-{source}.txids.txt",
    output:
        stats="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.stats.tsv",
        report=report(
            "results/{assembly}-{source}/{seqrepo}/report/mehari_db_check.txt",
            category="{assembly}-{source}",
            subcategory="{seqrepo}",
        ),
    log:
        "logs/{assembly}-{source}/mehari/{seqrepo}/check.log",
    conda:
        "../envs/datastuff.yaml"
    script:
        "../scripts/check_mehari_db.py"


rule generate_table_with_discards_to_review:
    input:
        discarded="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.report.jsonl",
    output:
        table="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.discards.tsv",
    log:
        "logs/{assembly}-{source}/mehari/{seqrepo}/discards_table.log",
    conda:
        "../envs/datastuff.yaml"
    shell:
        """(
        echo -e "kind\tidtype\tid\tgene_name\treason" > {output.table}
        jq -r 'select( .type == "Discard" and .value.source == "protobuf" ) | [.value.kind, .value.id.type, .value.id.value, .value.gene_name, .value.reason] | @tsv' {input.discarded} >> {output.table}
        ) >{log} 2>&1"""


rule datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/report.datavzrd.yaml"),
        mehari_check_db_stats="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.stats.tsv",
        fix_incorrect_cds="results/{assembly}-{source}/report/fix_incorrect_cds.tsv",
        refseq_id_to_ensembl_id="results/{assembly}-{source}/lookup/refseq_id_to_ensembl_id.tsv",
        discards_of_interest="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.discards.tsv",
    params:
        extra="",
    output:
        report(
            directory("results/{assembly}-{source}/datavzrd-report/{seqrepo}"),
            htmlindex="index.html",
            category="{assembly}-{source}",
            subcategory="{seqrepo}",
        ),
    log:
        "logs/datavzrd_report/{assembly}-{source}/{seqrepo}/check_mehari_db.log",
    wrapper:
        "v3.13.1/utils/datavzrd"
