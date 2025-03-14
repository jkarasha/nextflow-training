#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `INDEX` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

process QUANTIFICATION {
    tag "Salmon on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    //This seems a bit magical
    //It uses the mix operator to combine the two channels
    //and collect operator to collect all the files
    //Let's see what the mix operator does first
    //First let's check the quant channel:
    //quant_ch.view{ "Quantification: ${it}"}
    //Quantification: ~/dev/nextflow/nextflow-training/nf-training/work/6b/8d3ced2ce79b4c2206b881f0146903/liver
    //Quantification: ~/dev/nextflow/nextflow-training/nf-training/work/c6/50e44c4efeaa76b2a1dc7d531f45da/lung
    //Quantification: ~/dev/nextflow/nextflow-training/nf-training/work/10/b7c4dcbd6ce30505fc1cf5ae529b93/gut

    //Now let's check the fastqc channel:
    //fastqc_ch.view{ "FastQC: ${it}"}
    //FastQC: ~/dev/nextflow/nextflow-training/nf-training/work/d4/e15864edd9eaf97af795b21710ac9f/fastqc_gut_logs
    //FastQC: ~/dev/nextflow/nextflow-training/nf-training/work/f2/389e058c2ff9c1721f36745d2f703b/fastqc_liver_logs
    //FastQC: ~/dev/nextflow/nextflow-training/nf-training/work/16/a36163fb2562e5ce48adf4f1e2712f/fastqc_lung_logs

    //Now let's check the mix channel:
    //quant_ch.mix(fastqc_ch).view{ "Mix: ${it}"}
    //Mix: ~/dev/nextflow/nextflow-training/nf-training/work/16/a36163fb2562e5ce48adf4f1e2712f/fastqc_lung_logs
    //Mix: ~/dev/nextflow/nextflow-training/nf-training/work/d4/e15864edd9eaf97af795b21710ac9f/fastqc_gut_logs
    //Mix: ~/dev/nextflow/nextflow-training/nf-training/work/f2/389e058c2ff9c1721f36745d2f703b/fastqc_liver_logs
    //Mix: ~/dev/nextflow/nextflow-training/nf-training/work/10/b7c4dcbd6ce30505fc1cf5ae529b93/gut
    //Mix: ~/dev/nextflow/nextflow-training/nf-training/work/c6/50e44c4efeaa76b2a1dc7d531f45da/lung
    //Mix: ~/dev/nextflow/nextflow-training/nf-training/work/6b/8d3ced2ce79b4c2206b881f0146903/liver

    // Now lets's collect
    quant_ch.mix(fastqc_ch).collect().view{ "Collected: ${it}"}

    //collect() does the opposite of flatten(). It will take the 6 items and collect them into a single list
    //which explains why MULTIQC will run only once - a single task
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
}
