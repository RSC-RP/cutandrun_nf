include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main.nf'

//Samtools view to filter the mapped reads. 
workflow samtools_filter {
    take: 
    bam_bai_ch
    fasta

    main: 
    // Empty channel for the query name file used as a filter for now
    Channel.value( [] )
        .set{ qname }
    SAMTOOLS_VIEW(bam_bai_ch, fasta, qname)
    //Sort and index the picard marked dups/ filtered bams 
    SAMTOOLS_SORT(SAMTOOLS_VIEW.out.bam)
    // Index the picard marked dups/filtered bams
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    SAMTOOLS_SORT.out.bam
        .cross(SAMTOOLS_INDEX.out.bai){ row -> row[0].id } // join by the key name "id"
        .map { row -> [ row[0][0], row[0][1], row[1][1] ] }
        .set { bam_bai_ch }

    emit:
    bams_sorted     = SAMTOOLS_SORT.out.bam
    bam_bai_ch      = bam_bai_ch
}