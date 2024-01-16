include { FASTQC as FASTQC_TRIM } from './modules/nf-core/fastqc' //FASTQC_TRIM
include { TRIMGALORE } from './modules/nf-core/trimgalore'

workflow trimgalore {
    take:
    fastqs

    main:
    TRIMGALORE(meta_ch)
    FASTQC_TRIM(TRIMGALORE.out.reads)

    emit: 
    reads       =   TRIMGALORE.out.reads
    log         =   TRIMGALORE.out.log
    fastqc      =   FASTQC_TRIM.out.fastqc
}