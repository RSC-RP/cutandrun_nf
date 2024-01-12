include { MACS2_CALLPEAK } from '../../modules/nf-core/macs2/callpeak/main.nf'
include { KHMER_UNIQUEKMERS } from '../../modules/nf-core/khmer/uniquekmers/main.nf'

// Run MAC2 peak calling (subworkflow)
workflow macs2_peaks {
    take:
    bams
    fasta
    no_control
    effective_gsize
    calc_effective_gsize
    read_length

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
    if ( calc_effective_gsize ) {
        Channel.value(read_length)
            .set { khmer_size }
        KHMER_UNIQUEKMERS(fasta, khmer_size)
        KHMER_UNIQUEKMERS.out.gsize
            .set { gsize }
    } else {
        Channel.value(effective_gsize)
                .set { gsize }
    }
    // Run Macs2 peak calling
    MACS2_CALLPEAK(macs_ch, gsize)

    emit:
    macs2               =   MACS2_CALLPEAK.out.peak
    macs_ver            =   MACS2_CALLPEAK.out.versions
    khmer_ver           =   params.run_khmer ? KHMER_UNIQUEKMERS.out.versions : ''
}