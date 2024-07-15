rule initialize_seqrepo:
    input:
        ensembl="results/{assembly}-ensembl/reference/{assembly}-ensembl.fasta.gz",
        refseq="results/{assembly}-refseq/reference/{assembly}-refseq.fasta.gz",
        missing_fasta="results/{assembly}-{source}/fixes/missing.fasta",
    output:
        seqrepo_root=directory("results/{assembly}-{source}/seqrepo"),
        seqrepo_instance=directory("results/{assembly}-{source}/seqrepo/master"),
    params:
        ensembl_namespace=config["namespaces"]["ensembl"],
        refseq_namespace=config["namespaces"]["refseq"],
    log:
        "logs/{assembly}-{source}/seqrepo/initialize.log",
    conda:
        "../envs/seqrepo.yaml"
    shell:
        """
        seqrepo --root-directory {output.seqrepo_root} init --instance-name master

        # seqrepo load is too verbose and we cannot silence it
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root} \
            load --instance-name master --namespace {params.ensembl_namespace} \
            {input.ensembl} \
        | tail

        # seqrepo load is too verbose and we cannot silence it
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root} \
            load --instance-name master --namespace {params.refseq_namespace} \
            {input.refseq} \
        | tail

        # include missing records in both namespaces,
        # as seqrepo fetch somehow does not always seem to be working correctly (regarding the namespace)
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root} \
          load --instance-name master --namespace {params.refseq_namespace} {input.missing_fasta} | tail
        >{log} 2>&1 seqrepo --root-directory {output.seqrepo_root} \
          load --instance-name master --namespace {params.ensembl_namespace} {input.missing_fasta} | tail
        """


rule list_contigs:
    input:
        fasta="results/{assembly}-{source}/reference/{assembly}-{source}.fasta.gz",
    output:
        contigs="results/{assembly}-{source}/reference/{assembly}-{source}.contigs",
    log:
        "logs/{assembly}-{source}/seqrepo/list-contigs.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """ samtools dict -H {input.fasta} | cut -f2 | sed 's/^SN://' > {output.contigs} 2> {log}"""


rule detect_missing_sequences:
    input:
        unpack(cdot_input_mapping),
        ensembl="results/{assembly}-ensembl/reference/{assembly}-ensembl.contigs",
        refseq="results/{assembly}-refseq/reference/{assembly}-refseq.contigs",
    output:
        missing_txt=report(
            "results/{assembly}-{source}/mehari/seqrepo/missing.txt",
            category="{assembly}-{source}",
            subcategory="seqrepo",
        ),
        cdot_contigs_tmp=temp(
            "results/{assembly}-{source}/mehari/seqrepo/missing.txt.tmp"
        ),
        contigs_tmp=temp("results/{assembly}-{source}/mehari/seqrepo/contigs.txt.tmp"),
    params:
        cdot_files=lambda _, input: list(
            filter(
                lambda c: c is not None,
                (getattr(input, x, None) for x in ("cdot", "cdot_mt", "cdot_grafted")),
            )
        ),
    log:
        "logs/{assembly}-{source}/seqrepo/detect-missing.log",
    conda:
        "../envs/jq.yaml"
    shell:
        # or `.value.reason | contains("MissingSequence")`, but if there are more reasons a transcript has been discarded,
        # it makes no sense to fetch transcripts with missing sequence if they get discarded for a different reason anyway
        """(
        touch {output.cdot_contigs_tmp}
        for f in {params.cdot_files}; do
          jq -r '.transcripts[] | select(.id != null).id' <(pigz -dcf ${{f}}) >> {output.cdot_contigs_tmp}
        done;
        cat {output.cdot_contigs_tmp} | sort | uniq > {output.missing_txt}
        mv {output.missing_txt} {output.cdot_contigs_tmp}
        cat {input.ensembl} {input.refseq} | sort | uniq > {output.contigs_tmp}
        comm -23 {output.cdot_contigs_tmp} {output.contigs_tmp} > {output.missing_txt}
        )>{log} 2>&1"""


rule aggregate_missing_sequences:
    input:
        missing_accessions="results/{assembly}-{source}/mehari/seqrepo/missing.txt",
    output:
        missing_fasta="results/{assembly}-{source}/fixes/missing.fasta",
        missing_fasta_fai="results/{assembly}-{source}/fixes/missing.fasta.fai",
    params:
        missing_accessions=parse_missing_accessions,
    resources:
        ratelimit=1,
    log:
        "logs/{assembly}-{source}/mehari/seqrepo/fetch-missing.log",
    conda:
        "../envs/efetch.yaml"
    script:
        "../scripts/fetch_missing_sequences.py"
