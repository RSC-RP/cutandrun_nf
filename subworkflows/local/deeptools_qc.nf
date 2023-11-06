include { DEEPTOOLS_MULTIBIGWIGSUMMARY } from '../../modules/local/deeptools/multibigwigsummary.nf'
include { DEEPTOOLS_PLOTCORRELATION } from '../../modules/nf-core/deeptools/plotcorrelation/main.nf'
include { DEEPTOOLS_PLOTPCA } from '../../modules/nf-core/deeptools/plotpca/main.nf'
include { DEEPTOOLS_PLOTENRICHMENT as SEACR_PLOTENRICHMENT ; DEEPTOOLS_PLOTENRICHMENT as MACS2_PLOTENRICHMENT } from '../../modules/local/deeptools/plotenrichment.nf'

//Run deeptools bamcoverage to create signal/coverage file
workflow deeptools_qc {
    take: 
    bam_bai_ch
    bigwigs
    seacr_peaks
    run_macs2
    macs_peaks
    sample_sheet_name

    main: 
    // channel to collect all bigwigs in a single row
    bigwigs.map { row -> row[1] }
        .collect()
        .map { bw -> [ [ id:sample_sheet_name ], bw ] }
        .set { bw_ch }
    // channel to join bams and peak beds by meta.id per sample. 
    bam_bai_ch.cross(seacr_peaks){ row -> row[0].id } // join by the key name "id"
        .map { row -> [ row[0][0], [row[0][1]], [row[0][2]], [row[1][1]] ] }
        .set { bam_seacr_ch }

    // bam_peaks_ch.view{"the bam,bai, bed channel is $it"}
    // To compute the average values for all samples' coverage of the genome
    DEEPTOOLS_MULTIBIGWIGSUMMARY(bw_ch)

    // Use the genome-wide values to plot sample correlations and PCA
    channel.value([]).set{ method }
    channel.value([]).set{ plot_type }
    DEEPTOOLS_PLOTCORRELATION(DEEPTOOLS_MULTIBIGWIGSUMMARY.out.npz, method, plot_type)
    DEEPTOOLS_PLOTPCA(DEEPTOOLS_MULTIBIGWIGSUMMARY.out.npz)

    // calculate FRiP on all called peaks 
    SEACR_PLOTENRICHMENT(bam_seacr_ch)

    if ( run_macs2 ){
        bam_bai_ch.cross(macs_peaks){ row -> row[0].id } // join by the key name "id"
            .map { row -> [ row[0][0], [row[0][1]], [row[0][2]], [row[1][1]] ] }
            .set { bam_macs_ch }
        MACS2_PLOTENRICHMENT(bam_macs_ch)
    }

    emit:
    npz           = DEEPTOOLS_MULTIBIGWIGSUMMARY.out.npz
}