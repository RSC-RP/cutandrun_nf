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

    // Determine if IgG control bam will be used in the peak calling
    if ( no_control ) {
        bam_groups.target
            .map { meta, target_bam -> [ meta, target_bam, [] ] } // empty list for control bam 
            .set { macs_ch }
    } else {
        bam_groups.control
            .cross(bam_groups.targets){ meta -> meta[0].sample }
            .map { row -> [ row[1][0], row[1][1], row[0][1] ] } // [ meta, target_bam path, contro_bam path ]
            .set { macs_ch }
    }
    macs_ch.view { "the macs2 bam channel is $it"}
    // Use params.gsize for effective genome size or Run khmer to determine effective genome size
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
    khmer_ver           =   calc_effective_gsize ? KHMER_UNIQUEKMERS.out.versions : ''
}