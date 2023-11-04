include { DEEPTOOLS_MULTIBIGWIGSUMMARY } from '../../modules/local/deeptools/multibigwigsummary.nf'
include { DEEPTOOLS_PLOTCORRELATION } from '../../modules/nf-core/deeptools/plotcorrelation/main.nf'
include { DEEPTOOLS_PLOTPCA } from '../../modules/nf-core/deeptools/plotpca/main.nf'

//Run deeptools bamcoverage to create signal/coverage file
workflow deeptools_qc {
    take: 
    bigwigs
    sample_sheet_name

    main: 
    
    bigwigs.map { row -> row[1] }
        .collect()
        .map { bw -> [ [ id:sample_sheet_name ], bw ] }
        .set { bw_ch }
    DEEPTOOLS_MULTIBIGWIGSUMMARY(bw_ch)

    channel.value([]).set{ method }
    channel.value([]).set{ plot_type }
    DEEPTOOLS_PLOTCORRELATION(DEEPTOOLS_MULTIBIGWIGSUMMARY.out.npz, method, plot_type)
    DEEPTOOLS_PLOTPCA(DEEPTOOLS_MULTIBIGWIGSUMMARY.out.npz)

    emit:
    npz           = DEEPTOOLS_MULTIBIGWIGSUMMARY.out.npz
}