include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main.nf'
include { DEEPTOOLS_BAMCOVERAGE } from '../modules/nf-core/deeptools/bamcoverage/main.nf'

//Run deeptools bamcoverage to create signal/coverage file
workflow coverage_tracks {
    take: 
    bams
    fasta

    main: 
    //Sort and index the picard marked dups bams 
    SAMTOOLS_SORT(bams)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    //Create a channel for bedtools genomecov
    Channel.value( [] )
        .set { fai } //dummy channel for the fasta index file (fai)
    SAMTOOLS_SORT.out.bam
        .cross(SAMTOOLS_INDEX.out.bai){ meta -> meta[0].id } // join by the key name "id"
        .map { meta -> [ meta[0][0], meta[0][1], meta[1][1] ] }
        .set { coverage_ch }
    // calculate coverage track with Deeptools
    DEEPTOOLS_BAMCOVERAGE(coverage_ch, fasta, fai)
}