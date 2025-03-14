# Overview

## Transcriptome Index

DNA is transcribed into RNA, which is then translated into proteins. Different cells in your body have different parts of your DNA expressed. If you think of your DNA as a cookbook, the transcriptome is the set of recipes that are currently being used. The transcriptome is always changing, and different cells have different transcriptomes. Your transcriptome changes when you are sick vs when you are healthy. Your transcriptome change throughout the day.

A transcriptome index file is a reference data structure used in `RNA Sequencing` analysis to enable efficient alignment of RNA-seq reads to a set of transcript sequences.

They are created from reference transcriptome sequences(typically FASTA) and contain an optimized data structure that allows alignment tools to quickly locate where sequencing reads match in the transcriptome.

Some of the common tools used to create these include `bowtie2`, `hisat2`, `salmon`, and `kallisto`.

Transcriptome indexing is particularly important for RNA-seq because:
- It allows for efficient quantification of transcript abundance
- It enables detection of alternative splicing events
- It helps handle the complexity of mapping reads that might span exon-exon junctions

The index files themselves are typically binary and not meant to be read directly by humans, but they're essential intermediates in RNA-seq analysis pipelines.

## Channel Factories

`fromFilePairs`: Takes a glob pattern an input and returns a channel of tuples. The tuple has two items. First, is the read pair prefix, and the second is a list of paths to the read files.

For example, if the input is `data/ggal/gut_{1,2}.fq`, the output will be:

```console
[liver,
    [~/dev/nextflow/nextflow-training/nf-training/data/ggal/liver_1.fq,
     ~/dev/nextflow/nextflow-training/nf-training/data/ggal/liver_2.fq]
]
[gut,
    [~/dev/nextflow/nextflow-training/nf-training/data/ggal/gut_1.fq,
     ~/dev/nextflow/nextflow-training/nf-training/data/ggal/gut_2.fq]
]
[lung,
    [~/dev/nextflow/nextflow-training/nf-training/data/ggal/lung_1.fq,
     ~/dev/nextflow/nextflow-training/nf-training/data/ggal/lung_2.fq]
]
```

Channel factories also have options that can be used to modify their behavior. For example `checkIfExists` option can be used to check if the specified path contains file pairs. If the path does not contain file pairs, an error is thrown.

## Expression Quantification

Expression quantification is the process of measuring how active or "expressed" genes are in a biological sample. Determines how much each gene is being used by cells at a given time. When a gene is active or "expressed" your cell makes RNA copies of that gene. The more active the gene, the more RNA copies get made. Expression quantification is basically counting these RNA copies to see which genes are working hard and which ones are taking it easy or taking a break.


## Seqera Labs
