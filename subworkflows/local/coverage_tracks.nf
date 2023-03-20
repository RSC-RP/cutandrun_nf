include { SAMTOOLS_STATS } from '../../modules/nf-core/samtools/stats/main.nf'
include { DEEPTOOLS_BAMCOVERAGE } from '../../modules/nf-core/deeptools/bamcoverage/main.nf'

//Run deeptools bamcoverage to create signal/coverage file
workflow coverage_tracks {
    take: 
    bam_bai_ch
    fasta
    fai

    main: 
    // Calculate alignment QC stats
    SAMTOOLS_STATS(bam_bai_ch, fasta)
    // calculate coverage track with Deeptools
    DEEPTOOLS_BAMCOVERAGE(bam_bai_ch, fasta, fai)

    emit:
    stats           = SAMTOOLS_STATS.out.stats
    // bam_bai_ch      = bam_bai_ch
}