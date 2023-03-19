include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_STATS } from '../../modules/nf-core/samtools/stats/main.nf'
include { DEEPTOOLS_BAMCOVERAGE } from '../../modules/nf-core/deeptools/bamcoverage/main.nf'

//Run deeptools bamcoverage to create signal/coverage file
workflow coverage_tracks {
    take: 
    bams
    fasta
    fai

    main: 
    //Sort and index the picard marked dups bams 
    SAMTOOLS_SORT(bams)
    // Index the picard marked dups bams
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    SAMTOOLS_SORT.out.bam
        .cross(SAMTOOLS_INDEX.out.bai){ meta -> meta[0].id } // join by the key name "id"
        .map { meta -> [ meta[0][0], meta[0][1], meta[1][1] ] }
        .set { bam_bai_ch }
    // Calculate alignment QC stats
    SAMTOOLS_STATS(bam_bai_ch, fasta)
    // calculate coverage track with Deeptools
    DEEPTOOLS_BAMCOVERAGE(bam_bai_ch, fasta, fai)

    emit:
    stats = SAMTOOLS_STATS.out.stats
}