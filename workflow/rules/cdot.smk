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


rule fetch_incorrect_entries:
    params:
        transcripts=["NM_001137667.2", "NM_001137668.2", "NM_012115.4", "NM_001177639.3", "NM_001291281.3", "NM_001345921.3"]
    output:
        cdot="results/for-fix/{alias}/cdot.json.gz",
        xml="results/for-fix/{alias}/nuccore.xml",
    script:
        "../scripts/cdot_fetch_incorrect_entries.py"
