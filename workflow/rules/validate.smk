rule dump_mehari_db:
    input:
        tx_db="results/{alias}/mehari/{seqrepo}/txs.bin.zst",
    output:
        db_yaml="results/{alias}/mehari/{seqrepo}/txs.bin.zst.yaml.gz",
    # conda:
    #     "../envs/mehari.yaml"
    log:
        "logs/{alias}/mehari/{seqrepo}/dump.log",
    shell:
        """
        (
            mehari db dump --path-db {input.tx_db} | gzip -c > {output.db_yaml}
        ) >{log} 2>&1
        """


rule mehari_transcript_ids:
    input:
        db_yaml="results/{alias}/mehari/{seqrepo}/txs.bin.zst.yaml.gz",
    output:
        tx_ids="results/{alias}/mehari/{seqrepo}/txs.bin.zst.transcript_ids.txt",
    log:
        "logs/{alias}/mehari/{seqrepo}/transcript_ids.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """yq -r '.txDb.transcripts[].id' <(pigz -dc {input.db_yaml}) > {output.tx_ids} 2> {log}"""


rule cdot_tx_ids:
    input:
        unpack(cdot_input_mapping),
    output:
        tx_ids="results/{alias}/cdot/{alias}.txids.txt",
    log:
        "logs/{alias}/cdot/txids.log",
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
        kept="results/{alias}/mehari/{seqrepo}/txs.bin.zst.transcript_ids.txt",
        db_discarded="results/{alias}/mehari/{seqrepo}/txs.bin.zst.report.jsonl",
        cdot_tx_ids="results/{alias}/cdot/{alias}.txids.txt",
    output:
        stats="results/{alias}/mehari/{seqrepo}/txs.bin.zst.stats.tsv",
        report=report(
            "results/{alias}/{seqrepo}/report/mehari_db_check.txt",
            category="{alias}",
            subcategory="{seqrepo}",
        ),
    log:
        "logs/{alias}/mehari/{seqrepo}/check.log",
    conda:
        "../envs/datastuff.yaml"
    script:
        "../scripts/check_mehari_db.py"


rule datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/report.datavzrd.yaml"),
        mehari_check_db_stats="results/{alias}/mehari/{seqrepo}/txs.bin.zst.stats.tsv",
        fix_incorrect_entries="results/{alias}/report/fix_incorrect_entries.tsv",
        refseq_id_to_ensembl_id="results/{alias}/lookup/refseq_id_to_ensembl_id.tsv",
    params:
        extra="",
    output:
        report(
            directory("results/{alias}/datavzrd-report/{seqrepo}/"),
            htmlindex="index.html",
            category="{alias}",
            subcategory="{seqrepo}",
        ),
    log:
        "logs/datavzrd_report/{alias}/{seqrepo}/check_mehari_db.log",
    wrapper:
        "v3.13.1/utils/datavzrd"
