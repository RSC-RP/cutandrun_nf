include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main.nf'
include { PICARD_MARKDUPLICATES; PICARD_MARKDUPLICATES as PICARD_RMDUPLICATES } from '../../modules/nf-core/picard/markduplicates/main.nf'
include { SAMTOOLS_STATS } from '../../modules/nf-core/samtools/stats/main.nf'
include { DEEPTOOLS_BAMCOVERAGE } from '../../modules/nf-core/deeptools/bamcoverage/main.nf'

//Run deeptools bamcoverage to create signal/coverage file
workflow markdups_bigwigs {
    take: 
    bams
    fasta
    fai

    main: 
    //Sort and index the picard marked dups bams 
    SAMTOOLS_SORT(bams)
    PICARD_MARKDUPLICATES(SAMTOOLS_SORT.out.bam, fasta , fai)
    if ( params.remove_dups ){
        PICARD_RMDUPLICATES(SAMTOOLS_SORT.out.bam)
        PICARD_RMDUPLICATES.out.bam
            .set { bams_sorted }
    } else {
        PICARD_MARKDUPLICATES.out.bam
            .set { bams_sorted }
    }
    // Index the picard marked dups bams
    SAMTOOLS_INDEX(bams_sorted)
    bams_sorted
        .cross(SAMTOOLS_INDEX.out.bai){ meta -> meta[0].id } // join by the key name "id"
        .map { meta -> [ meta[0][0], meta[0][1], meta[1][1] ] }
        .set { bam_bai_ch }
    // Calculate alignment QC stats
    SAMTOOLS_STATS(bam_bai_ch, fasta)
    // calculate coverage track with Deeptools
    DEEPTOOLS_BAMCOVERAGE(bam_bai_ch, fasta, fai)

    emit:
    bams = SAMTOOLS_SORT.out.bam
    metrics = PICARD_MARKDUPLICATES.out.metrics
    stats = SAMTOOLS_STATS.out.stats
}