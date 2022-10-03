nextflow.enable.dsl = 2

//Include Modules
include { TRIMGALORE } from './modules/nf-core/modules/trimgalore/main.nf'
include { BOWTIE2_ALIGN; BOWTIE2_ALIGN as SPIKEIN_ALIGN } from './modules/nf-core/modules/bowtie2/align/main.nf'
include { BAMTOBEDGRAPH } from './modules/local/bedtools/main.nf'
include { SEACR_CALLPEAK } from './modules/nf-core/modules/seacr/callpeak/main.nf'
include { PICARD_MARKDUPLICATES; PICARD_MARKDUPLICATES as PICARD_RMDUPLICATES } from './modules/nf-core/modules/picard/markduplicates/main.nf'
include { SAMTOOLS_SORT; SAMTOOLS_SORT as SAMTOOLS_NSORT } from './modules/nf-core/modules/samtools/sort/main.nf'
include { DEEPTOOLS_BAMCOVERAGE } from './modules/nf-core/modules/deeptools/bamcoverage/main.nf'

//Include subworkflows
include { bowtie2_index; bowtie2_index as bowtie2_index_spike } from './subworkflows/bowtie_index.nf'
include { macs2_peaks } from './subworkflows/macs2_peaks.nf'

//Define stdout message for the command line use
idx_or_fasta = (params.index == '' ? params.fasta : params.index)
log.info """\
         C U T & R U N-  P I P E L I N E
         ===================================
         Project           : $workflow.projectDir
         Project workDir   : $workflow.workDir
         Container Engine  : $workflow.containerEngine
         Samples           : ${params.sample_sheet}
         Genome            : ${idx_or_fasta}
         """
         .stripIndent()

//run the workflow for alignment, to bedgraph, to peak calling for Cut&Run data
workflow call_peaks {
        //Empty channel to collect the versions of software used 
        Channel.empty()
            .set { versions }
        //Stage the fasta file(s)
        Channel.fromPath(file(params.fasta, checkIfExists: true))
            .collect()
            .set { fasta }
        //Optionally, create the index from a fasta file
        if ( params.build_index ) {
            bowtie2_index(fasta)
            bowtie2_index.out.index
                .collect() //collect converts this to a value channel and to be used multiple times
                .set { index }
            versions = versions.concat(bowtie2_index.out.versions)
        } else {
            //Stage the genome index directory
            Channel.fromPath(file(params.index, checkIfExists: true))
                .collect()
                .set { index }
        }
        //Optionally, create the spike-in index from a fasta file
        if ( params.build_spike_index ) {
            //Stage the fasta files
            Channel.fromPath(file(params.spike_fasta, checkIfExists: true))
                .set { spike_fasta }
            bowtie2_index_spike(spike_fasta)
            bowtie2_index_spike.out.index
                .collect()
                .set { spike_index }
            versions = versions.concat(bowtie2_index_spike.out.versions)
        } else {
            //Stage the genome index directory
            Channel.fromPath(file(params.spike_index, checkIfExists: true))
                .collect() //collect converts this to a value channel and used multiple times
                .set { spike_index }
        }
        //Create the input channel which contains the SAMPLE_ID, whether its single-end, and the file paths for the fastqs. 
        Channel.fromPath(file(params.sample_sheet, checkIfExists: true))
            .splitCsv(header: true, sep: ',')
            .map { meta -> [ [ "id":meta["sample_id"], "single_end":meta["single_end"].toBoolean(), "group":meta["target_or_control"], "sample":meta["sample"] ], //meta
                             [ file(meta["read1"], checkIfExists: true), file(meta["read2"], checkIfExists: true) ] //reads
                           ] }
            .set { meta_ch }
        //Stage the file for bedgraph generations
        Channel.fromPath(file(params.chrom_sizes, checkIfExists: true))
            .collect()
            .set { chrom_sizes }
        //Add fastqc module here 
        //Adapter and Quality trimming of the fastq files 
        TRIMGALORE(meta_ch)
        //Perform the alignement
        spike_in = false
        //NOTE: should have bowtie2 save unaligned reads to a seperate file. 
        BOWTIE2_ALIGN(TRIMGALORE.out.reads, index, spike_in,
                      params.save_unaligned)
        if ( params.spike_norm ){
            spike_in = true
            SPIKEIN_ALIGN(TRIMGALORE.out.reads, spike_index, spike_in,
                         params.save_unaligned)
            //Channel containing tuples of the spikeIn seqdepth and scaling factor constant "C"
            Channel.value(params.scale_factor_constant)
                .set { C }
            SPIKEIN_ALIGN.out.seq_depth
                .combine( C )
                .map { val -> [ "seq_depth":val[0].toInteger() , "constant":val[1] ] }
                .map { val -> if (val.seq_depth > 0 ) {
                                    val.constant.div(val.seq_depth) 
                                } else {
                                    val.constant.div(1)
                                } }
                .set { scale_factor }
        } else {
            //If not using the spikeIn normalization, then just need empty list
            Channel.value( [] )
                .set { scale_factor }
        }

        //Add samtools index module here
        //Add picard markDuplicates modules here
        PICARD_MARKDUPLICATES(BOWTIE2_ALIGN.out.bam)
        if ( params.remove_dups ){
            PICARD_RMDUPLICATES(BOWTIE2_ALIGN.out.bam)
            PICARD_RMDUPLICATES.out.bam
                .set { bams }
        } else {
            PICARD_MARKDUPLICATES.out.bam
                .set { bams }
        }
        //Sort bam files by read names
        SAMTOOLS_NSORT(bams)
        //Add Samtools stats module here for QC
        //Add [optional] samtools quality score filtering here 
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
        //SEACR peak calling
        SEACR_CALLPEAK(seacr_ch, params.threshold)
        // MACS2 peak calling, Optional
        if ( params.run_macs2 ){
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
            //Run MAC2 peak calling
            macs2_peaks(macs_ch)
        }
        // calculate coverage track with Deeptools
        SAMTOOLS_SORT(bams)
        Channel.value( [[]] )
            .set { bai } //dummy channel for the bam index file (bai)
        Channel.value( [] )
            .set { fai } //dummy channel for the fasta index file (fai)
        SAMTOOLS_SORT.out.bam
            .combine( bai )
            .set { coverage_ch }
        DEEPTOOLS_BAMCOVERAGE(coverage_ch, fasta, fai)

        //Add Deeptools module to calculate FRIP here
        //Add multiQC module here 
        // versions.concat(TRIMGALORE.out.versions, 
        //                 BOWTIE2_ALIGN.out.versions,
        //                 BAMTOBEDGRAPH.out.versions, 
        //                 SEACR_CALLPEAK.out.versions)
        //         .collect()
        //         .collectFile(name: 'versions.txt', newLine: true)
}

//End with a message to print to standard out on workflow completion. 
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

