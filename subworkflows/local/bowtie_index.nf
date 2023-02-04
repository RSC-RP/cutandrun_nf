include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main.nf'

//Generate the index file (subworkflow)
workflow bowtie2_index {
    take:
    fasta

    main:
    //execute the BOWTIE2 genome index process
    BOWTIE2_BUILD(fasta)

    emit:
    index = BOWTIE2_BUILD.out.index
    versions = BOWTIE2_BUILD.out.versions
}