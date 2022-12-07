include { BAMTOBEDGRAPH } from '../modules/local/bedtools/main.nf'
include { SEACR_CALLPEAK } from '../modules/nf-core/modules/seacr/callpeak/main.nf'
include { SAMTOOLS_SORT as SAMTOOLS_NSORT } from '../modules/nf-core/modules/samtools/sort/main.nf'

// Run SEACR peak caller
workflow seacr_peaks {
    take:
    bams
    chrom_sizes
    scale_factor

    main:
    //Sort bam files by read names for bam to bedpe conversion
    SAMTOOLS_NSORT(bams)
    // Convert the bam files to bed format with bedtools 
    Channel.value(params.spike_norm)
        .set { spike_norm }
    BAMTOBEDGRAPH(SAMTOOLS_NSORT.out.bam, chrom_sizes, spike_norm, scale_factor)
    //Split the control and the target channels. should be a custom groovy function really - need to figure this out.
    BAMTOBEDGRAPH.out.bedgraph
        .branch { 
            control: it[0].group =~ /control/
            targets: it[0].group =~ /target/
        }
        .set { bedgraphs }
    //Define SEACR formatted input channels
    if ( params.threshold > 0 ){
        //Channel containing an empty/dummy value for the control file
        bedgraphs.targets
            .combine( Channel.value( [[]] ) )
            .set { seacr_ch }
    } else {
        //Channel containing the targets and the control bedgraphs 
        bedgraphs.control
            .cross(bedgraphs.targets){ meta -> meta[0].sample } // join by the key name "sample"
            .map { meta -> [ meta[1][0], meta[1][1], meta[0][1] ] }
            .set { seacr_ch }
    }
    //Convert the numeric threshold to a value channel
    Channel.value(params.threshold)
        .set { threshold }
    //SEACR peak calling
    SEACR_CALLPEAK(seacr_ch, threshold, spike_norm)
}

