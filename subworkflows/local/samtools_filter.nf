include { SAMTOOLS_STATS } from '../../modules/nf-core/samtools/stats/main.nf'

//Samtools view to filter the mapped reads. 
workflow samtools_filter {
    take: 
    bam_bai_ch
    fasta

    main: 
    // Calculate alignment QC stats
    // SAMTOOLS_STATS(bam_bai_ch, fasta)
    // calculate coverage track with Deeptools
    SAMTOOLS_VIEW(bam_bai_ch)

    emit:
    bam = SAMTOOLS_VIEW.out.bam
}