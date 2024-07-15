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
            mehari db dump --path-db {input.tx_db} | pigz -c > {output.db_yaml}
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
        """yq -r '.txDb.transcripts[] | select(.filtered == false) | .id' <(pigz -dc {input.db_yaml}) > {output.tx_ids} 2> {log}"""


rule mehari_hgnc_ids:
    input:
        db_yaml="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.yaml.gz",
    output:
        hgnc_ids="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.hgnc_ids.txt",
    log:
        "logs/{assembly}-{source}/mehari/{seqrepo}/hgnc_ids.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """yq -r '.txDb.geneToTx[] | select(.filtered == false) | .geneId' <(pigz -dc {input.db_yaml}) > {output.hgnc_ids} 2> {log}"""


rule cdot_tx_ids:
    input:
        unpack(cdot_input_mapping),
    output:
        transcript_ids="results/{assembly}-{source}/cdot/{assembly}-{source}.transcript_ids.txt",
    log:
        "logs/{assembly}-{source}/cdot/transcript_ids.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """
        (
        touch {output.transcript_ids}
        for f in {input}; do
            pigz -dcf ${{f}} | jq -r '.transcripts | keys[]' >> {output.transcript_ids}
        done
        ) >{log} 2>&1
        """


rule cdot_hgnc_ids:
    input:
        unpack(cdot_input_mapping),
    output:
        hgnc_ids="results/{assembly}-{source}/cdot/{assembly}-{source}.hgnc_ids.txt",
        tmp=temp(
            "results/{assembly}-{source}/cdot/{assembly}-{source}.hgnc_ids.txt.tmp"
        ),
    log:
        "logs/{assembly}-{source}/cdot/hgnc_ids.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """
        (
        touch {output.hgnc_ids}; touch {output.tmp}
        for f in {input}; do
            pigz -dcf ${{f}} | jq -r '.genes[].hgnc | select(. != null)' >> {output.tmp}
        done
        cat {output.tmp} | sort -n | uniq > {output.hgnc_ids}
        ) >{log} 2>&1
        """


rule cdot_gene_symbols:
    input:
        unpack(cdot_input_mapping),
    output:
        hgnc_ids="results/{assembly}-{source}/cdot/{assembly}-{source}.gene_symbols.txt",
        tmp=temp(
            "results/{assembly}-{source}/cdot/{assembly}-{source}.gene_symbols.txt.tmp"
        ),
    log:
        "logs/{assembly}-{source}/cdot/gene_symbols.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """
        (
        touch {output.hgnc_ids}; touch {output.tmp}
        for f in {input}; do
            pigz -dcf ${{f}} | jq -r '.genes[].gene_symbol | select(. != null)' >> {output.tmp}
        done
        cat {output.tmp} | sort | uniq > {output.hgnc_ids}
        ) >{log} 2>&1
        """


rule check_mehari_db:
    input:
        unpack(cdot_input_mapping),
        tx_db="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst",
        tx_db_report="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.report.jsonl",
        hgnc="results/hgnc/hgnc_complete_set.json",
        genes_to_disease="results/human-phenotype-ontology/genes_to_disease_with_hgnc_id.tsv",
    output:
        # stats="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.stats.tsv",
        report=report(
            "results/{assembly}-{source}/mehari/{seqrepo}/report/mehari_db_check.txt",
            category="{assembly}-{source}",
            subcategory="{seqrepo}",
        ),
    params:
        known_issues=get_known_issues,
        cdot=get_mehari_cdot_param_string,
    log:
        "logs/{assembly}-{source}/mehari/{seqrepo}/check.log",
    # conda:
    #     "../envs/mehari.yaml"
    shell:
        """(
        mehari db check \
        --path-db {input.tx_db} \
        --path-hgnc-json {input.hgnc} \
        --path-disease-gene-tsv {input.genes_to_disease} \
        {params.cdot} \
        --output {output.report} \
        ) >{log} 2>&1"""


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
        echo -e "idtype\tid\tgene_name\treason\ttags" > {output.table}
        jq -r 'select( .type == "Discard" ) | [.value.id.type, .value.id.value, .value.gene_name, .value.reason, (.value.tags // [] | join(","))] | @tsv' {input.discarded} >> {output.table}
        ) >{log} 2>&1"""


rule datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/report.datavzrd.yaml"),
        mehari_check_db_stats="results/{assembly}-{source}/mehari/{seqrepo}/txs.bin.zst.stats.tsv",
        fix_incorrect_cds="results/{assembly}-{source}/cdot/fix_incorrect_cds.tsv",
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
