nextflow.enable.dsl = 2

// Include Modules
include { FASTQC; FASTQC as FASTQC_TRIM } from './modules/nf-core/fastqc'
include { TRIMGALORE } from './modules/nf-core/trimgalore'
include { MULTIQC } from './modules/nf-core/multiqc'
include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx'
include { BOWTIE2_ALIGN; BOWTIE2_ALIGN as SPIKEIN_ALIGN } from './modules/nf-core/bowtie2/align'
include { PICARD_MARKDUPLICATES; PICARD_MARKDUPLICATES as PICARD_RMDUPLICATES } from './modules/nf-core/picard/markduplicates/main.nf'
include { SAMTOOLS_STATS } from './modules/nf-core/samtools/stats/main.nf'

// Include subworkflows
include { samtools_filter } from './subworkflows/local/samtools_filter.nf'
include { coverage_tracks } from './subworkflows/local/coverage_tracks.nf'
include { deeptools_qc } from './subworkflows/local/deeptools_qc.nf'
include { seacr_peaks } from './subworkflows/local/seacr_peaks.nf'
include { macs2_peaks } from './subworkflows/local/macs2_peaks.nf'
include { bowtie2_index; bowtie2_index as bowtie2_index_spike } from './subworkflows/local/bowtie2_index.nf'

// Define stdout message for the command line use
idx_or_fasta = (params.index == '' ? params.fasta : params.index)
log.info """\
         C U T & R U N-  P I P E L I N E
         ===================================
         Project           : $workflow.projectDir
         Project workDir   : $workflow.workDir
         Container Engine  : $workflow.containerEngine
         Results           : ${params.outdir}
         Samples           : ${params.sample_sheet}
         Genome            : ${idx_or_fasta}
         """
         .stripIndent()

// Run the workflow for alignment, to bedgraph, to peak calling for Cut&Run data
workflow align_call_peaks {
        //Empty channel to collect the versions of software used 
        Channel.empty()
            .set { versions }
        //Stage the and index fasta file(s)
        Channel.fromPath(file(params.fasta, checkIfExists: true))
            .map { fasta -> [ [], fasta ] } // fasta files now need meta info
            .collect()
            .set { fasta }
        Channel.fromPath(file(params.spike_fasta, checkIfExists: true))
                .map { fasta -> [ [], fasta ] } // fasta files now need meta info
                .collect()
                .set { spike_fasta }
        SAMTOOLS_FAIDX(fasta)
        SAMTOOLS_FAIDX.out.fai
            .collect()
            .set { fai }
        //Optionally, create the bowtie2 index from a fasta file
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
            bowtie2_index_spike(spike_fasta)
            bowtie2_index_spike.out.index
                .collect()
                .set { spike_index }
            versions = versions.concat(bowtie2_index_spike.out.versions)
        } else if ( params.spike_norm ) {
            //Stage the genome index directory
            Channel.fromPath(file(params.spike_index, checkIfExists: true))
                .collect() //collect converts this to a value channel and used multiple times
                .set { spike_index }
        }
        //Stage the file for bedgraph generations
        Channel.fromPath(file(params.chrom_sizes, checkIfExists: true))
            .collect()
            .set { chrom_sizes }
        //Stage the multiqc configurations
        Channel.fromPath(file(params.multiqc_config,checkIfExists: true))
            .collect()
            .set { multiqc_config }
        //Create the input channel which contains the SAMPLE_ID, whether its single-end, and the file paths for the fastqs. 
        Channel.fromPath(file(params.sample_sheet, checkIfExists: true))
            .splitCsv(header: true, sep: ',')
            .map { meta -> [ [ "id":meta["sample_id"], "single_end":meta["single_end"].toBoolean(), "group":meta["target_or_control"], "sample":meta["sample"] ], //meta
                             [ file(meta["read1"], checkIfExists: true), file(meta["read2"], checkIfExists: true) ] //reads
                           ] }
            .set { meta_ch }
        // define static variables 
        def sample_sheet_name = file(params.sample_sheet, checkIfExists: true).simpleName
        def run_macs2 = params.run_macs2
        def effective_gsize = params.gsize
        def no_control = params.no_control_macs2
        def calc_effective_gsize = params.calc_effective_gsize
        def read_length = params.read_length

        //fastqc of raw sequence
        FASTQC(meta_ch)
        //Adapter and Quality trimming of the fastq files 
        TRIMGALORE(meta_ch)
        FASTQC_TRIM(TRIMGALORE.out.reads)

        // Perform the alignement
        spike_in = false
        BOWTIE2_ALIGN(TRIMGALORE.out.reads, index, spike_in,
                      params.save_unaligned)
        if ( params.spike_norm ){
            spike_in = true
            SPIKEIN_ALIGN(TRIMGALORE.out.reads, spike_index, spike_in,
                         params.save_unaligned)
            // Channel for BAMs aligned to target genome, spike-in seq depth, and calculated scale_factor
            Channel.value(params.scale_factor_constant)
                .set { C }
            SPIKEIN_ALIGN.out.seq_depth
                .combine( C )
                .map { meta, seq_depth, C -> 
                        depth = seq_depth.toInteger()
                        if ( depth > 0 ) {
                            sf = C.div(depth) 
                        } else {
                            sf = C.div(1)
                        }
                        meta.spike_seq_depth = seq_depth
                        meta.scale_factor = sf
                        return [ meta, [] ]
                    }
                .cross(BOWTIE2_ALIGN.out.bam) { meta -> meta[0].id }
                .map { sf_meta, bam -> [ sf_meta[0], bam[1] ]}
                .set { bams_ch }

        } else {
            BOWTIE2_ALIGN.out.bam
                .set { bams_ch }
            //If not using the spikeIn normalization, then just need empty lists
            Channel.value( [] )
                .set { scale_factor }
        }

        // Run picard markduplicates, optionally remove duplicates
        // note: Bowtie2 modules automatically sorts the input BAMs to picard
        PICARD_MARKDUPLICATES(bams_ch, fasta, fai)
        //create bam and bai channel for samtool stats
        PICARD_MARKDUPLICATES.out.bam
            .cross(PICARD_MARKDUPLICATES.out.bai){ row -> row[0].id } // join by the key name "id"
            .map { row -> [ row[0][0], row[0][1], row[1][1] ] }
            .set { mkdup_ch }
        // Calculate alignment QC stats on the markedDup bams 
        SAMTOOLS_STATS(mkdup_ch, fasta)

        //Optionally, remove duplicates before peak calling
        if ( params.remove_dups ){
            PICARD_RMDUPLICATES(BOWTIE2_ALIGN.out.bam, fasta, fai)
            //create bam and bai channel
            PICARD_RMDUPLICATES.out.bam
                .cross(PICARD_RMDUPLICATES.out.bai){ row -> row[0].id } // join by the key name "id"
                .map { row -> [ row[0][0], row[0][1], row[1][1] ] }
                .set { bam_bai_ch }
            //create channel with sorted bams (no bai)
            PICARD_RMDUPLICATES.out.bam
                .set { bams_sorted }
        } else {
            //rename channel to standardize inputs to peak callers
            mkdup_ch
                .set { bam_bai_ch }
            //create channel with sorted bams (no bai)
            PICARD_MARKDUPLICATES.out.bam
                .set { bams_sorted }
        }

        //Optionally: samtools quality score and/or alignment flag filtering 
        if ( params.filter_bam ){
            samtools_filter(bam_bai_ch, fasta)
            // replace bam and bai channel
            samtools_filter.out.bam_bai_ch
                .set { bam_bai_ch }
            // replace channel with sorted bams
            samtools_filter.out.bams_sorted
                .set { bams_sorted }
        }

        // Create coverage (bigwig or bedgraph) files for IGV/UCSC and look at genomic enrichment stats
        coverage_tracks(bam_bai_ch, fasta, fai, sample_sheet_name)

        // SEACR peak calling 
        seacr_peaks(bams_sorted, chrom_sizes)

        // MACS2 peak calling, Optional
        if ( run_macs2 ){
            //Run MAC2 peak calling
            macs2_peaks(bams_sorted, fasta, no_control, effective_gsize, calc_effective_gsize, read_length)
            macs2_peaks.out.macs2 
                .set { macs_peaks }
        } else {
            // set a dummy channel if MACS is not run
            Channel.value([])
                .set { macs_peaks }
        }

        // PCA, Correlation, and FRiP for samples
        deeptools_qc(bam_bai_ch, 
                     coverage_tracks.out.bigwig, 
                     seacr_peaks.out.seacr, 
                     run_macs2,
                     macs_peaks, 
                     sample_sheet_name)

        //MultiQC to collect QC results
        if ( params.spike_norm ) { 
            SPIKEIN_ALIGN.out.log
                .set { spike_log }
        } else {
            Channel.value([])
                .set { spike_log }
        }

        //Create channel for all the QC metrics to be included in MultiQC
        FASTQC.out.fastqc
            .concat(TRIMGALORE.out.log)
            .concat(FASTQC_TRIM.out.fastqc)
            .concat(BOWTIE2_ALIGN.out.log)
            .concat(spike_log)
            .concat(PICARD_MARKDUPLICATES.out.metrics)
            .concat(SAMTOOLS_STATS.out.stats)
            .concat(coverage_tracks.out.metrics)
            .concat(coverage_tracks.out.matrix)
            .concat(deeptools_qc.out.corr)
            .concat(deeptools_qc.out.pca)
            .map { row -> row[1] }
            .collect()
            .set { multiqc_ch }

        if (params.extra_multiqc_config){
            Channel.fromPath(file(params.extra_multiqc_config, checkIfExists: true))
                .collect()
                .set { extra_multiqc_config }
        }else{
            Channel.value([])
                .set{ extra_multiqc_config }
        }
        if (params.multiqc_logo){
            Channel.fromPath(file(params.multiqc_logo, checkIfExists:true))
                .set { multiqc_logo }
        }else{
            Channel.value([])
                .set { multiqc_logo }
        }
        // MultiQC process
        MULTIQC(multiqc_ch, 
                multiqc_config,
                extra_multiqc_config, 
                multiqc_logo, 
                sample_sheet_name)
}


workflow call_peaks {
    //Create the input channel for bams which contains the SAMPLE_ID, whether its single-end, and the file paths. 
    Channel.fromPath(file(params.sample_sheet, checkIfExists: true))
        .splitCsv(header: true, sep: ',')
        .map { row -> [ [ "id":row["sample_id"], "single_end":row["single_end"].toBoolean(), "group":row["target_or_control"], "sample":row["sample"] ], //meta
                        [ file(row["bam"], checkIfExists: true) ] //reads
                    ] }
        .set { bams }
    //Stage the file for bedgraph generations
    Channel.fromPath(file(params.chrom_sizes, checkIfExists: true))
        .collect()
        .set { chrom_sizes }
    //I need a scale factor per bam if spikein norm has been run previously. 
    if ( params.spike_norm ) {
        //Create a channel for the scale factor for each BAM file in the sample sheet
        Channel.fromPath(file(params.sample_sheet, checkIfExists:true))
            .map { row -> row["scale_factor"] }
            .set { scale_factor }
    } else {
        Channel.value( [] )
            .set { scale_factor }
    }
    // SEACR peak calling
    seacr_peaks(bams, chrom_sizes, scale_factor)
    // MACS2 peak calling, Optional
    if ( params.run_macs2 ) {
        //Run MAC2 peak calling
        macs2_peaks(bams)
    }
}

workflow bowtie2_index_only {
    //Stage the fasta file(s)
    Channel.fromPath(file(params.fasta, checkIfExists: true))
        .collect()
        .set { fasta }

    // Create the index from a fasta file
    bowtie2_index(fasta)
    bowtie2_index.out.index
        .collect() //collect converts this to a value channel and to be used multiple times
        .set { index }

    //Optionally create the spike-in index 
    if ( params.build_spike_index ) {
        //Stage the fasta files
        Channel.fromPath(file(params.spike_fasta, checkIfExists: true))
            .collect()
            .set { spike_fasta }
        bowtie2_index_spike(spike_fasta)
        bowtie2_index_spike.out.index
            .collect()
            .set { spike_index }
    }
}

//End with a message to print to standard out on workflow completion. 
// workflow.onComplete {
//     println "Pipeline completed at: $workflow.complete"
//     if ( $workflow.success ){
//         println "Execution status: OK "
//     } else {
//         println "Execution status: ${workflow.errorMessage}"
//     }
// }
