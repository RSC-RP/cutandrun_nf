include { DEEPTOOLS_BAMCOVERAGE } from '../../modules/nf-core/deeptools/bamcoverage/main.nf'
include { DEEPTOOLS_PLOTFINGERPRINT } from '../../modules/nf-core/deeptools/plotfingerprint/main.nf'

//Run deeptools bamcoverage to create signal/coverage file
workflow coverage_tracks {
    take: 
    bam_bai_ch
    fasta
    fai
    sample_sheet_name

    main: 
    // channel for all samples in a single plot
    // concatenate with the individual samples 
    bam_bai_ch
        .map { meta, bam, bai ->  [  "bam", bam ,
                                     "bai", bai ] }
        .groupTuple( by: [0,2] )
        .map { row -> [ [ id:sample_sheet_name ], row[1], row[3] ] }
        .concat( bam_bai_ch )
        .set { fingerprint_ch }
    
    // calculate coverage track with Deeptools
    DEEPTOOLS_BAMCOVERAGE(bam_bai_ch, fasta, fai)

    // Investigate enrichment of genomic regions
    DEEPTOOLS_PLOTFINGERPRINT(fingerprint_ch)

    emit:
    bigwig              =   DEEPTOOLS_BAMCOVERAGE.out.bigwig
    metrics             =   DEEPTOOLS_PLOTFINGERPRINT.out.metrics
    matrix              =   DEEPTOOLS_PLOTFINGERPRINT.out.matrix
}