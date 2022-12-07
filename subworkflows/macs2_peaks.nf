include { MACS2_CALLPEAK } from '../modules/nf-core/modules/macs2/callpeak/main.nf'
include { KHMER_UNIQUEKMERS } from '../modules/nf-core/modules/khmer/uniquekmers/main.nf'

// Run MAC2 peak calling (subworkflow)
workflow macs2_peaks {
    take:
    bams

    main:
    //Separate the bam files by target or control antibody
    bams.branch { 
            control: it[0].group =~ /control/
            targets: it[0].group =~ /target/
        }
        .set { bam_groups }
    bam_groups.control
        .cross(bam_groups.targets){ meta -> meta[0].sample }
        .map { meta -> [ meta[1][0], meta[1][1], meta[0][1] ] }
        .set { macs_ch }
    //Either run khmer to determine effective genome size for macs2, or use a value provided as params.gsize
    if ( params.run_khmer ) {
        Channel.fromPath(file(params.fasta, checkIfExists: true))
            .set { fasta }
        Channel.value(params.kmer_size)
            .set { khmer_size }
        KHMER_UNIQUEKMERS(fasta, khmer_size)
        KHMER_UNIQUEKMERS.out.kmers
                .set { gsize }
    } else {
        Channel.value(params.gsize)
                .set { gsize }
    }
    // Run Macs2 peak calling
    MACS2_CALLPEAK(macs_ch, gsize)

    emit:
    macs_ver = MACS2_CALLPEAK.out.versions
    khmer_ver = params.run_khmer ? KHMER_UNIQUEKMERS.out.versions : ''
}