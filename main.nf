nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from './modules/nf-core/modules/bowtie2/build/main.nf'
include { TRIMGALORE } from './modules/nf-core/modules/trimgalore/main.nf'
include { BOWTIE2_ALIGN } from './modules/nf-core/modules/bowtie2/align/main.nf'
include { BAMTOBEDGRAPH } from './modules/local/bedtools/main.nf'
include { SEACR_CALLPEAK } from './modules/nf-core/modules/seacr/callpeak/main.nf'
include { MACS2_CALLPEAK } from './modules/nf-core/modules/macs2/callpeak/main.nf'

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
//QUESTION: need a switch of MACS vs SEACR peak caller to be used? or do both all the time?
workflow call_peaks {
        //Empty channel to collect the versions of software used 
        Channel.empty()
            .set { versions }
        //if there is no filepath to the index provided, create the index from a fasta file
        if ( params.build_index ) {
            bowtie2_index()
            bowtie2_index.out.index
                .set { index }
            versions = versions.concat(bowtie2_index.out.versions)
        } else {
            //Stage the genome index directory
            Channel.fromPath(file(params.index, checkIfExists: true))
                .collect() //collect converts this to a value channel and used multiple times
                .set { index }
        }
        //Create the input channel which contains the SAMPLE_ID, whether its single-end, and the file paths for the fastqs. 
        Channel.fromPath(file(params.sample_sheet, checkIfExists: true))
            .splitCsv(header: true, sep: ',')
            .map { meta -> [ [ "id":meta["sample_id"], "single_end":meta["single_end"].toBoolean(), "group":meta["target_or_control"], "sample":meta["sample"] ], //meta
                             [ file(meta["read1"], checkIfExists: true), file(meta["read2"], checkIfExists: true) ] //reads
                           ] }
            .set { meta_ch }
        //Stage the file for bedgraph generations
        Channel.fromPath(file(params.genome_file, checkIfExists: true))
            .collect()
            .set { genome_file }
        //Adapter and Quality trimming of the fastq files 
        TRIMGALORE(meta_ch)
        //Perform the alignement 
        BOWTIE2_ALIGN(TRIMGALORE.out.reads, index,
                      params.save_unaligned, params.sort_bam)
        //Conver the bam files to bed format with bedtools 
        BAMTOBEDGRAPH(BOWTIE2_ALIGN.out.bam, genome_file)
        //Define the control and the target channels. should be a custom groovy function really - need to figure this out.
        BAMTOBEDGRAPH.out.bedgraph
            .filter( ~/^.+control.+/ )
            .set { controls }
        BAMTOBEDGRAPH.out.bedgraph
            .filter( ~/^.+target.+/ )
            .set { targets }
        controls
            .cross(targets){ meta -> meta[0].sample }
            .map { meta -> [ meta[1][0], meta[1][1], meta[0][1] ] }
            .set { seacr_ch }
        
        //SEACR peak calling
        SEACR_CALLPEAK(seacr_ch, params.threshold)
        //MACS2 peak calling 
        BOWTIE2_ALIGN.out.bam
            .branch { 
               control: it =~ /^.+control.+/
               targets: it =~ /^.target.+/
            }
            .set { conditions }
        conditions.control
            .cross(conditions.targets){ meta -> meta[0].sample }
            .map { meta -> [ meta[1][0], meta[1][1], meta[0][1] ] }
            .set { macs_ch }
        // MACS2_CALLPEAK(macs_ch, params.macs2_gsize)

        // versions.concat(TRIMGALORE.out.versions, 
        //                 BOWTIE2_ALIGN.out.versions,
        //                 BAMTOBEDGRAPH.out.versions)
        //         .collect()
        //         .collectFile(name: 'versions.txt', newLine: true)
}

//Generate the index file 
workflow bowtie2_index {
    main:        
    //Stage the fasta files
    Channel.fromPath(params.fasta)
        .ifEmpty { error "No fasta files found ${params.fasta}. Required to build genome index with bowtie2" }
        .set { fasta }
    //execute the BOWTIE2 genome index process
    BOWTIE2_BUILD(fasta)

    emit:
    index = BOWTIE2_BUILD.out.index
    versions = BOWTIE2_BUILD.out.versions
}

//End with a message to print to standard out on workflow completion. 
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

