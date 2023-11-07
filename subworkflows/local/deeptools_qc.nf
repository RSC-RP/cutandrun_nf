include { DEEPTOOLS_MULTIBIGWIGSUMMARY } from '../../modules/local/deeptools/multibigwigsummary.nf'
include { DEEPTOOLS_PLOTCORRELATION } from '../../modules/nf-core/deeptools/plotcorrelation/main.nf'
include { DEEPTOOLS_PLOTPCA } from '../../modules/nf-core/deeptools/plotpca/main.nf'
include { DEEPTOOLS_PLOTENRICHMENT as SEACR_PLOTENRICHMENT ; DEEPTOOLS_PLOTENRICHMENT as MACS2_PLOTENRICHMENT } from '../../modules/local/deeptools/plotenrichment.nf'
include { MACSPEAKSTOBED } from '../../modules/local/macspeakstobed.nf'

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

    // To compute the average values for all samples' coverage of the genome
    DEEPTOOLS_MULTIBIGWIGSUMMARY(bw_ch)

    // Use the genome-wide values to plot sample correlations and PCA
    channel.value([]).set{ method }
    channel.value([]).set{ plot_type }
    DEEPTOOLS_PLOTCORRELATION(DEEPTOOLS_MULTIBIGWIGSUMMARY.out.npz, method, plot_type)
    DEEPTOOLS_PLOTPCA(DEEPTOOLS_MULTIBIGWIGSUMMARY.out.npz)

    // calculate FRiP on called peaks per sample
    SEACR_PLOTENRICHMENT(bam_seacr_ch)
    if ( run_macs2 ){
        // convert narrow/broad/gapped peaks file to standard 3 column bed file
        MACSPEAKSTOBED(macs_peaks)
        // join the bams, bai, and peaks by meta.id
        bam_bai_ch.cross(MACSPEAKSTOBED.out.bed){ row -> row[0].id } // join by the key name "id"
            .map { row -> [ row[0][0], row[0][1], row[0][2], row[1][1] ] }
            .set { bam_macs_ch }
        // generate enrichment/FRiP stats per sample
        MACS2_PLOTENRICHMENT(bam_macs_ch)
    }

    emit:
    npz             =   DEEPTOOLS_MULTIBIGWIGSUMMARY.out.npz
    corr            =   DEEPTOOLS_PLOTCORRELATION.out.matrix
    pca             =   DEEPTOOLS_PLOTPCA.out.tab
}