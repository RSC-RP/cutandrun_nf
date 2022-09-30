include { MACS2_CALLPEAK } from '../modules/nf-core/modules/macs2/callpeak/main.nf'
include { KHMER_UNIQUEKMERS } from '../modules/nf-core/modules/khmer/uniquekmers/main.nf'

// Run MAC2 peak calling (subworkflow)
workflow macs2_peaks {
    take:
    macs_ch

    main:
    //Either run khmer to determine effective genome size for macs2, or use a value provided as params.macs2_gsize
    if ( params.run_khmer ) {
        Channel.fromPath(file(params.fasta, checkIfExists: true))
            .set { fasta }
        Channel.value(params.kmer_size)
            .set { khmer_size }
        KHMER_UNIQUEKMERS(fasta, khmer_size)
        KHMER_UNIQUEKMERS.out.kmers
                .set { macs2_gsize }
    } else {
        Channel.value(params.macs2_gsize)
                .set { macs2_gsize }
    }
    // Run Macs2 peak calling
    MACS2_CALLPEAK(macs_ch, macs2_gsize)

    emit:
    macs_ver = MACS2_CALLPEAK.out.versions
    khmer_ver = params.run_khmer ? KHMER_UNIQUEKMERS.out.versions : ''
}