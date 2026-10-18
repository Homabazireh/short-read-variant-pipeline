nextflow.enable.dsl = 2

process PREPARE_BENCHMARK {
    tag "${meta.key}"
    container params.containers.bcftools
    publishDir "${params.outdir}/callsets/${meta.key}", mode: 'copy'
    input:
    tuple val(meta), path(query, stageAs: 'query.input.vcf.gz'), path(truth, stageAs: 'truth.input.vcf.gz'), path(bed, stageAs: 'evaluation.bed')
    output:
    tuple val(meta), path('query.vcf.gz'), path('query.vcf.gz.tbi'), path('truth.vcf.gz'), path('truth.vcf.gz.tbi'), path('evaluation.bed')
    path 'eligibility.tsv'
    script:
    def select = meta.variant_type == 'sv' ? "--include 'N_ALT=1 && TYPE=\"indel\" && (STRLEN(REF)=1 || STRLEN(ALT)=1) && ABS(ILEN)>=50'" : ''
    """
    bcftools view ${select} query.input.vcf.gz -Oz -o query.vcf.gz
    bcftools view ${select} truth.input.vcf.gz -Oz -o truth.vcf.gz
    bcftools index --tbi query.vcf.gz
    bcftools index --tbi truth.vcf.gz
    printf 'role\tinput_records\teligible_records\n' > eligibility.tsv
    printf 'query\t%s\t%s\n' "\$(bcftools view -H query.input.vcf.gz | wc -l)" "\$(bcftools view -H query.vcf.gz | wc -l)" >> eligibility.tsv
    printf 'truth\t%s\t%s\n' "\$(bcftools view -H truth.input.vcf.gz | wc -l)" "\$(bcftools view -H truth.vcf.gz | wc -l)" >> eligibility.tsv
    """
}

process HAPPY {
    tag "${meta.key}"
    container params.containers.happy
    publishDir "${params.outdir}/callsets/${meta.key}", mode: 'copy'
    input:
    tuple val(meta), path(query, stageAs: 'query.vcf.gz'), path(query_tbi, stageAs: 'query.vcf.gz.tbi'), path(truth, stageAs: 'truth.vcf.gz'), path(truth_tbi, stageAs: 'truth.vcf.gz.tbi'), path(bed, stageAs: 'evaluation.bed')
    tuple path(fasta, stageAs: 'reference.fa'), path(fai, stageAs: 'reference.fa.fai')
    output:
    path 'happy.*'
    script:
    """
    mkdir -p tmp
    export TMPDIR="\$PWD/tmp"
    export RTG_MEM=${Math.max(1, (task.memory.toGiga() * 0.75).intValue())}g
    /opt/hap.py/bin/hap.py truth.vcf.gz query.vcf.gz --engine vcfeval --adjust-conf-regions --threads ${task.cpus} \
        --reference reference.fa -f evaluation.bed -o happy
    """
}

process TRUVARI {
    tag "${meta.key}"
    container params.containers.truvari
    publishDir "${params.outdir}/callsets/${meta.key}", mode: 'copy'
    input:
    tuple val(meta), path(query, stageAs: 'query.vcf.gz'), path(query_tbi, stageAs: 'query.vcf.gz.tbi'), path(truth, stageAs: 'truth.vcf.gz'), path(truth_tbi, stageAs: 'truth.vcf.gz.tbi'), path(bed, stageAs: 'evaluation.bed')
    tuple path(fasta, stageAs: 'reference.fa'), path(fai, stageAs: 'reference.fa.fai')
    output:
    path 'truvari'
    script:
    """
    truvari bench --base truth.vcf.gz --comp query.vcf.gz --output truvari \
        --reference reference.fa --includebed evaluation.bed \
        --bSample ${meta.truth_sample} --cSample ${meta.sample} \
        --refdist ${meta.matching.refdist} --pctseq ${meta.matching.pctseq} \
        --pctsize ${meta.matching.pctsize} --pctovl ${meta.matching.pctovl} \
        --sizemin ${meta.matching.sizemin} --sizefilt ${meta.matching.sizefilt} \
        --sizemax ${meta.matching.sizemax} --pick ${meta.matching.pick} \
        --passonly --no-ref a --extend 0
    """
}

workflow {
    jobs = Channel.fromList(params.jobs).map { row -> tuple(row, file(row.vcf), file(row.truth_vcf), file(row.evaluation_bed)) }
    reference = Channel.value(tuple(file(params.reference_fasta), file(params.reference_fai)))
    PREPARE_BENCHMARK(jobs)
    prepared = PREPARE_BENCHMARK.out[0].branch { meta, query, qi, truth, ti, bed ->
        small: meta.variant_type == 'small'
        sv: meta.variant_type == 'sv'
    }
    HAPPY(prepared.small, reference)
    TRUVARI(prepared.sv, reference)
}
