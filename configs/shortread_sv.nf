#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Sarek 3.10.0 disables TIDDIT assembly by passing an empty BWA index.
// Keep the assembly-capable SV branch explicit, isolated, and independently resumable.
process CLASSIC_BWA_INDEX {
    container 'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path 'bwa', emit: index
    path 'bwa.version.txt', emit: version

    script:
    """
    mkdir bwa
    bwa index -p bwa/reference.fa '${fasta}'
    bwa 2>&1 | sed -n 's/^Version: //p' > bwa.version.txt
    """
}

process TIDDIT_ASSEMBLY {
    tag "$sample"
    container 'community.wave.seqera.io/library/tiddit:3.9.5--3fb6c287f34e6ab0'
    publishDir { "${params.outdir}/raw/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path fasta, stageAs: 'reference.fa'
    path fai, stageAs: 'reference.fa.fai'
    path bwa

    output:
    tuple val(sample), path("${sample}.tiddit.raw.vcf"), emit: vcf
    path "${sample}.tiddit.raw.ploidies.tab", emit: ploidies
    path "${sample}.tiddit.version.txt", emit: version

    script:
    """
    for extension in amb ann bwt pac sa; do
        ln -s bwa/reference.fa.\${extension} reference.fa.\${extension}
    done
    tiddit --sv --threads ${task.cpus} --bam '${bam}' --ref reference.fa -o '${sample}.tiddit.raw'
    tiddit --help > '${sample}.tiddit.version.txt' 2>&1
    """
}

process RESTRICT_SV_REGIONS {
    tag "$sample"
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'
    publishDir { "${params.outdir}/variants/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(vcf)
    path regions

    output:
    tuple val(sample), path("${sample}.tiddit.vcf.gz"), path("${sample}.tiddit.vcf.gz.tbi"), emit: variants
    path "${sample}.bcftools.version.txt", emit: version

    script:
    """
    bcftools sort '${vcf}' -Oz -o sorted.vcf.gz
    bcftools index -t sorted.vcf.gz
    bcftools view -R '${regions}' --regions-overlap 1 sorted.vcf.gz -Oz -o '${sample}.tiddit.vcf.gz'
    bcftools index -t '${sample}.tiddit.vcf.gz'
    bcftools --version > '${sample}.bcftools.version.txt'
    """
}

workflow {
    if (!params.input || !params.fasta || !params.fai || !params.regions || !params.outdir) {
        error 'Required parameters: input, fasta, fai, regions, outdir'
    }
    alignments = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam, checkIfExists: true), file(row.bai, checkIfExists: true)) }
    CLASSIC_BWA_INDEX(file(params.fasta, checkIfExists: true))
    TIDDIT_ASSEMBLY(alignments, file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true), CLASSIC_BWA_INDEX.out.index)
    RESTRICT_SV_REGIONS(TIDDIT_ASSEMBLY.out.vcf, file(params.regions, checkIfExists: true))
    RESTRICT_SV_REGIONS.out.variants
        .map { sample, vcf, tbi -> "${sample},tiddit,${params.outdir}/variants/${sample}/${vcf.name},sv\n" }
        .collectFile(name: 'calls.csv', storeDir: params.outdir, sort: true,
                     seed: 'sample,label,vcf,variant_type\n', newLine: false)
}
