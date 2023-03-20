include { DEEPTOOLS_BAMCOVERAGE } from '../../modules/nf-core/deeptools/bamcoverage/main.nf'

//Run deeptools bamcoverage to create signal/coverage file
workflow coverage_tracks {
    take: 
    bam_bai_ch
    fasta
    fai

    main: 
    // calculate coverage track with Deeptools
    DEEPTOOLS_BAMCOVERAGE(bam_bai_ch, fasta, fai)

    emit:
    bigwig           = DEEPTOOLS_BAMCOVERAGE.out.bigwig
}