rule dump_mehari_db:
    input:
        tx_db="results/mehari/{alias}/{seqrepo}/txs.bin.zst",
    output:
        db_yaml="results/mehari/{alias}/{seqrepo}/txs.bin.zst.yaml.gz",
    # conda:
    #     "../envs/mehari.yaml"
    log:
        "logs/mehari/{alias}/{seqrepo}/dump.log",
    shell:
        """
        mehari db dump --path-db {input.tx_db} | gzip -c > {output.db_yaml} 2> {log}
        """


rule mehari_transcript_ids:
    input:
        db_yaml="results/mehari/{alias}/{seqrepo}/txs.bin.zst.yaml.gz",
    output:
        tx_ids="results/mehari/{alias}/{seqrepo}/txs.bin.zst.transcript_ids.txt",
    log:
        "logs/mehari/{alias}/{seqrepo}/transcript_ids.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """yq -r '.txDb.transcripts[].id' <(pigz -dc {input.db_yaml}) > {output.tx_ids} 2> {log}"""


rule check_mehari_db:
    input:
        kept="results/mehari/{alias}/{seqrepo}/txs.bin.zst.transcript_ids.txt",
        db_discarded="results/mehari/{alias}/{seqrepo}/txs.bin.zst.report.jsonl",
        cdot="results/transcripts/cdot/{alias}.hgnc.json.gz",
    output:
        stats=report(
            "results/mehari/{alias}/{seqrepo}/txs.bin.zst.stats.tsv",
            category="{alias}",
            subcategory="{seqrepo}",
        ),
        report=report(
            "results/report/{alias}/{seqrepo}/mehari_db_check.txt",
            category="{alias}",
            subcategory="{seqrepo}",
        ),
    conda:
        "../envs/datastuff.yaml"
    script:
        "../scripts/check_mehari_db.py"
