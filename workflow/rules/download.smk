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
    shell:
        """
        (
        for i in {{00..03}};
         do wget -qO- "https://github.com/varfish-org/clinvar-data-jsonl/releases/download/clinvar-weekly-20250908/clinvar-data-jsonl-20250908+0.18.5.tar.gz.$i";
        done
        ) | pigz -dc | tar -f- -xO "clinvar-data-jsonl-20250908+0.18.5/clinvar-full-release.jsonl.gz" > {output.clinvar}
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
        (
          echo -e "vcv\tscv\tsymbol\ttxAcc";
          pigz -dc {input.clinvar} \
          | jq -r 'select(.classifiedRecord.clinicalAssertions != null) | .accession as $vcv | .classifiedRecord.clinicalAssertions[] | .clinvarAccession.accession as $scv | select(.simpleAllele.genes != null and (.simpleAllele.genes | length) > 0) | (.simpleAllele.genes[0].symbol // "UNKNOWN") as $gene | .simpleAllele.attributes[]? | select(.attribute.type == "HGVS" and ((.attribute.base.value | startswith("NM_")) or (.attribute.base.value | startswith("NR_")))) | "\($vcv)\t\($scv)\t\($gene)\t\(.attribute.base.value)"' | cut -d : -f 1 | cut -d . -f 1
        ) | sed -e "s/SCV004027544.NM_177438/SCV004027544\tDICER1\tNM_177438/g" > {output.clinvar_info}

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

        (
          echo -e "hgncId";
          jq -r '[.classifiedRecord.simpleAllele.genes[]? | .hgncId?] | map(select(. != null and . != "")) | unique[]' <(pigz -dc {input.clinvar}) \
        ) > {output.hgnc_ids}

        mlr --itsv --otsv count -g hgncId then sort -nr count {output.hgnc_ids} > {output.hgnc_id_counts}
        """
