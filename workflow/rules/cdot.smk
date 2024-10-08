rule extract_tags:
    input:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.json.gz",
    output:
        tags="results/{assembly}-{source}/cdot/{assembly}-{source}.tags.tsv",
    log:
        "logs/{assembly}-{source}/cdot/extract_tags.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/extract_tags.py"


rule fetch_incorrect_entries:
    params:
        transcripts=transcripts_to_fix_with_nuccore,
    output:
        xml="results/{assembly}-{source}/cdot/nuccore.xml.gz",
    log:
        "logs/{assembly}-{source}/cdot/fetch_incorrect_entries.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/fetch_incorrect_entries.py"


rule fix_incorrect_cds:
    input:
        xml="results/{assembly}-{source}/cdot/nuccore.xml.gz",
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.hgnc.json.gz",
    output:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.hgnc.cds.json.gz",
        report="results/{assembly}-{source}/cdot/fix_incorrect_cds.tsv",
    conda:
        "../envs/datastuff.yaml"
    log:
        "logs/{assembly}-{source}/cdot/fix_incorrect_cds.log",
    script:
        "../scripts/cdot/fix_incorrect_cds.py"


rule update_from_hgnc_complete_set:
    input:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.json.gz",
        hgnc="results/hgnc/hgnc_complete_set.json",
    output:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.hgnc.json.gz",
        report=report(
            "results/report/{assembly}-{source}/cdot_hgnc_update.tsv",
            category="{assembly}-{source}",
        ),
    log:
        "logs/{assembly}-{source}/cdot/update_with_hgnc_complete_set.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cdot/update_from_hgnc.py"


rule extract_chrMT:
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


# RefSeq or MANE select transcripts
rule find_select_transcripts:
    input:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.hgnc.json.gz",
    output:
        accessions="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.select.txt",
    log:
        "logs/{assembly}-{source}/report/find_select_transcripts.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """(
        pigz -dc {input.cdot} | jq -r '.transcripts[] | select(.genome_builds[].tag) | select (. != null) | select( .genome_builds[].tag | contains("Select")).id' | sort > {output.accessions}
        ) >{log} 2>&1"""


rule find_partial_transcripts:
    input:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.hgnc.json.gz",
    output:
        accessions="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.partial.txt",
    log:
        "logs/{assembly}-{source}/report/find_partial_transcripts.log",
    conda:
        "../envs/jq.yaml"
    shell:
        """(
        pigz -dc {input.cdot} | jq -r '.transcripts[] | select(.partial == 1 ).id' | sort > {output.accessions}
        ) >{log} 2>&1"""


rule determine_partial_but_select_transcripts:
    input:
        select="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.select.txt",
        partial="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.partial.txt",
    output:
        accessions="results/{assembly}-{source}/cdot/{assembly}-{source}.partial_but_select.txt",
    conda:
        "../envs/base.yaml"
    log:
        "logs/{assembly}-{source}/report/determine_partial_but_select_transcripts.log",
    shell:
        """(comm -12 {input.select} {input.partial} > {output.accessions}) >{log} 2>&1"""


rule lookup_ensembl_ids_for_refseq_ids:
    input:
        accessions="results/{assembly}-{source}/cdot/{assembly}-{source}.partial_but_select.txt",
    params:
        additional_accessions=transcripts_to_lookup_ensembl_ids_for,
    output:
        tsv="results/{assembly}-{source}/lookup/refseq_id_to_ensembl_id.tsv",
    log:
        "logs/{assembly}-{source}/lookup/lookup_ensembl_ids_for_refseq_ids.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/lookup_ensembl_ids_for_refseq_ids.py"


rule cdot_graft_ensembl_ids_for_certain_refseq_ids:
    input:
        cdot=ensembl_cdot,
        lookup="results/{assembly}-{source}/lookup/refseq_id_to_ensembl_id.tsv",
    output:
        cdot="results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.grafted.json.gz",
    conda:
        "../envs/datastuff.yaml"
    log:
        "logs/{assembly}-{source}/cdot/graft_ensembl_ids.log",
    script:
        "../scripts/cdot/graft_ensembl_ids.py"
