rule get_ensembl_sequence:
    output:
        "results/{assembly}-ensembl/reference/{assembly}-ensembl.{datatype}.fasta",
    params:
        species=get_ensembl_sequence_param("species"),
        build=get_ensembl_sequence_param("build"),
        release=get_ensembl_sequence_param("release"),
        datatype=lambda wildcards: wildcards.datatype,
    log:
        "logs/{assembly}-ensembl/get_ensembl_sequence.{datatype}.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v5.0.2/bio/reference/ensembl-sequence"


rule merge_ensembl_sequence:
    input:
        "results/{assembly}-ensembl/reference/{assembly}-ensembl.cdna.fasta",
        "results/{assembly}-ensembl/reference/{assembly}-ensembl.ncrna.fasta",
    output:
        "results/{assembly}-ensembl/reference/{assembly}-ensembl.fasta.gz",
    log:
        "logs/{assembly}-ensembl/merge_sequence.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        (
        cat {input} | bgzip -c > {output}
        ) >{log} 2>&1
        """


rule get_refseq_installed_files:
    output:
        "results/{assembly}-refseq/reference/files.installed",
    params:
        species=get_refseq_sequence_param("species"),
        species_name=get_refseq_sequence_param("species_name"),
    log:
        "logs/{assembly}-refseq/get_refseq_installed_files.log",
    conda:
        "../envs/base.yaml"
    shell:
        """(
        curl --silent "https://ftp.ncbi.nih.gov/refseq/{params.species}/mRNA_Prot/{params.species_name}.files.installed" > {output}
        ) >{log} 2>&1"""


rule get_refseq_sequence:
    input:
        refseq_installed="results/{assembly}-refseq/reference/files.installed",
    output:
        refseq_seq="results/{assembly}-refseq/reference/{assembly}-refseq.fasta.gz",
    params:
        species=get_refseq_sequence_param("species"),
        species_name=get_refseq_sequence_param("species_name"),
    cache: "omit-software"
    log:
        "logs/{assembly}-refseq/get_refseq_sequence.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/download_refseq.py"


rule get_cdot_data:
    output:
        "results/{assembly}-{source}/cdot/{assembly}-{source}.cdot.json.gz",
    params:
        url=get_cdot_download_url,
    conda:
        "../envs/base.yaml"
    cache: "omit-software"
    log:
        "logs/{assembly}-{source}/get_cdot_transcripts.log",
    shell:
        """(wget --quiet -O - {params.url} | gzip -dcf | jq '.' | gzip -c > {output}) 2> {log}"""


rule get_hgnc_complete_set:
    output:
        "results/hgnc/hgnc_complete_set.json",
    params:
        url=get_hgnc_complete_set_download_url,
    conda:
        "../envs/base.yaml"
    cache: "omit-software"
    log:
        "logs/hgnc/get_hgnc_complete_set.log",
    shell:
        """(wget -O -  {params.url} | jq '.' > {output}) 2> {log}"""


rule get_genes_to_disease:
    output:
        "results/human-phenotype-ontology/genes_to_disease.tsv",
    params:
        url=get_genes_to_disease_download_url(),
    conda:
        "../envs/base.yaml"
    cache: "omit-software"
    log:
        "logs/human-phenotype-ontology/get_genes_to_disease.log",
    shell:
        """wget --quiet -O {output} {params.url} 2> {log}"""


rule add_hgnc_id_to_genes_to_disease:
    input:
        genes_to_disease="results/human-phenotype-ontology/genes_to_disease.tsv",
        hgnc_complete_set="results/hgnc/hgnc_complete_set.json",
    output:
        genes_to_disease="results/human-phenotype-ontology/genes_to_disease_with_hgnc_id.tsv",
    conda:
        "../envs/datastuff.yaml"
    log:
        "logs/human-phenotype-ontology/add_hgnc_id_to_genes_to_disease.log",
    script:
        "../scripts/human_phenotype_ontology/add_hgnc_id_to_genes_to_disease.py"


rule get_clinvar_affected_transcripts:
    output:
        clinvar="results/clinvar/clinvar.jsonl.gz",
    conda:
        "../envs/base.yaml"
    params:
        release=config["clinvar-data-jsonl"]["release"],
        date=config["clinvar-data-jsonl"]["release"].split("+", 1)[0],
    shell:
        """
        (
        for i in {{00..04}};
         do wget -qO- "https://github.com/varfish-org/clinvar-data-jsonl/releases/download/clinvar-weekly-{params.date}/clinvar-data-jsonl-{params.release}.tar.gz.$i";
        done
        ) | pigz -dc | tar -f- -xO "clinvar-data-jsonl-{params.release}/clinvar-full-release.jsonl.gz" > {output.clinvar}
        """


rule clinvar_tx_accs:
    input:
        clinvar="results/clinvar/clinvar.jsonl.gz",
    output:
        clinvar_info="results/clinvar/clinvar-vcv-scv-symbol-txAcc.tsv",
        tx_acc_count="results/clinvar/clinvar-txAcc-count.tsv",
    conda:
        "../envs/base.yaml"
    shell:
        """
        expression='select(.classifiedRecord?.clinicalAssertions? != null)
          | .accession as $vcv
          | .classifiedRecord.clinicalAssertions[]?
          | .clinvarAccession?.accession as $scv
          | (.simpleAllele? // {{}}) as $sa
          | select( ($sa.genes? // []) | length > 0)
          | ($sa.genes?[0]?.symbol? // "UNKNOWN") as $raw_gene
          # FIX: this record has `NM_177438.3:c.5527+26A>G` as its _gene symbol_
          | (if $scv == "SCV004027544" then "DICER1" else $raw_gene end) as $gene
          | $sa.attributes?[]?
          | select(.attribute?.type == "HGVS" and
                   ((.attribute?.base?.value? // "") | (startswith("NM_") or startswith("NR_"))))
          | (.attribute?.base?.value? | split(":")[0] | split(".")[0]) as $tx
          | "\($vcv)\t\($scv)\t\($gene)\t\($tx)"'
        (
          echo -e "vcv\tscv\tsymbol\ttxAcc";
          pigz -dc {input.clinvar} | jaq -r "$expression"
        ) > {output.clinvar_info}

        mlr --itsv --otsv count -g txAcc then sort -nr count {output.clinvar_info} > {output.tx_acc_count}
        """


rule clinvar_hgnc_id_counts:
    input:
        clinvar="results/clinvar/clinvar.jsonl.gz",
    output:
        hgnc_ids=temp("results/clinvar/clinvar-hgnc-ids.tsv"),
        hgnc_id_counts="results/clinvar/clinvar-hgnc-id-counts.tsv",
    conda:
        "../envs/base.yaml"
    shell:
        """
        expression='[ .classifiedRecord?.simpleAllele?.genes?[]? | .hgncId? ] | map(select(. != null and . != ""))'
        (
          echo -e "hgncId";
          pigz -dc {input.clinvar} | jaq -r "$expression" | sort -u
        ) > {output.hgnc_ids}

        mlr --itsv --otsv count -g hgncId then sort -nr count {output.hgnc_ids} > {output.hgnc_id_counts}
        """
