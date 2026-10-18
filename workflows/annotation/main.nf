nextflow.enable.dsl = 2

process PREPARE {
    tag "${meta.key}"
    container params.containers.bcftools
    publishDir "${params.outdir}/callsets/${meta.key}", mode: 'copy'
    input:
    tuple val(meta), path(vcf, stageAs: 'input.vcf.gz')
    tuple path(fasta, stageAs: 'reference.fa'), path(fai, stageAs: 'reference.fa.fai')
    output:
    tuple val(meta), path("${meta.key}.prepared.vcf.gz"), path("${meta.key}.prepared.vcf.gz.tbi")
    script:
    def normalize = meta.variant_type == 'small' ? 'norm --check-ref e -m -any -f reference.fa' : 'view'
    """
    bcftools ${normalize} input.vcf.gz -Ou | bcftools sort -Oz -o ${meta.key}.prepared.vcf.gz
    bcftools index --tbi ${meta.key}.prepared.vcf.gz
    """
}

process VEP {
    tag "${meta.key}"
    container params.containers.vep
    publishDir "${params.outdir}/callsets/${meta.key}", mode: 'copy', pattern: '*.html'
    input:
    tuple val(meta), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi')
    tuple path(fasta, stageAs: 'reference.fa'), path(fai, stageAs: 'reference.fa.fai')
    path cache, stageAs: 'cache'
    tuple path(clinvar, stageAs: 'clinvar.vcf.gz'), path(clinvar_index, stageAs: 'clinvar.vcf.gz.tbi')
    tuple path(gnomad, stageAs: 'gnomad.vcf.gz'), path(gnomad_index, stageAs: 'gnomad.vcf.gz.tbi')
    output:
    tuple val(meta), path("${meta.key}.vep.vcf")
    path '*.html', optional: true
    script:
    def custom = []
    if (meta.variant_type == 'small' && params.annotation.clinvar) {
        custom << "--custom file=clinvar.vcf.gz,short_name=CLINVAR,format=vcf,type=exact,coords=0,fields=${params.annotation.clinvar.fields.join('%')}"
    }
    if (meta.variant_type == 'small' && params.annotation.gnomad) {
        custom << "--custom file=gnomad.vcf.gz,short_name=GNOMAD,format=vcf,type=exact,coords=0,fields=${params.annotation.gnomad.fields.join('%')}"
    }
    if (meta.records == 0) {
        """
        gzip -dc input.vcf.gz > ${meta.key}.vep.vcf
        """
    } else {
        """
        vep --offline --cache --dir_cache cache --cache_version 116 --assembly GRCh38 \
            --species homo_sapiens --fasta reference.fa --format vcf --vcf \
            --input_file input.vcf.gz --output_file ${meta.key}.vep.vcf \
            --stats_file ${meta.key}.vep.html --fork ${task.cpus} --force_overwrite \
            --symbol --hgvs --mane --canonical --biotype --numbers --variant_class \
            --protein --flag_pick --allele_number ${custom.join(' ')}
        """
    }
}

process COMPRESS {
    tag "${meta.key}"
    container params.containers.bcftools
    publishDir "${params.outdir}/callsets/${meta.key}", mode: 'copy'
    input:
    tuple val(meta), path(vcf, stageAs: 'annotated.vcf')
    output:
    tuple val(meta), path("${meta.key}.annotated.vcf.gz"), path("${meta.key}.annotated.vcf.gz.tbi")
    script:
    """
    bcftools view annotated.vcf -Oz -o ${meta.key}.annotated.vcf.gz
    bcftools index --tbi ${meta.key}.annotated.vcf.gz
    """
}

workflow {
    jobs = Channel.fromList(params.jobs).map { row -> tuple(row, file(row.vcf)) }
    reference = Channel.value(tuple(file(params.reference_fasta), file(params.reference_fai)))
    cache = Channel.value(file(params.annotation.vep_cache))
    clinvar = Channel.value(params.annotation.clinvar ? tuple(file(params.annotation.clinvar.path), file(params.annotation.clinvar.index)) : tuple([], []))
    gnomad = Channel.value(params.annotation.gnomad ? tuple(file(params.annotation.gnomad.path), file(params.annotation.gnomad.index)) : tuple([], []))
    PREPARE(jobs, reference)
    VEP(PREPARE.out, reference, cache, clinvar, gnomad)
    COMPRESS(VEP.out[0])
}
