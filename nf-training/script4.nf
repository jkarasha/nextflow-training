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

    tag "Running Salmon on $sample_id"
    publishDir params.outdir, mode: 'copy'

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

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    read_pairs_ch.view { it }

    index_ch = INDEX(params.transcriptome_file)
    // Notice how the quantification process is defined with two inputs
    // The first is the path to the "index" and the second is the tuple of read pairs, id included
    // For each read pair, the process will be executed once.args
    // If pass in a glob pattern, it will be expanded into a list of paths

    // nextflow run script4.nf -resume --reads 'data/ggal/*_{1,2}.fq' will get all read pairs in ggal folder
    // The glob pattern will be expanded into a three tuples, which will trigger three tasks
    // Really important concept in nextflow. If the input channel has {n} items,
    // it will cause the process to be spawn {n} number of tasks. The tasks will run in parallel.
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)

    // Print the quantification results
    quant_ch.println()
}
