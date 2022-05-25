nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from './modules/nf-core/modules/bowtie2/build/main.nf'
include { TRIMGALORE } from './modules/nf-core/modules/trimgalore/main.nf'
include { BOWTIE2_ALIGN } from './modules/nf-core/modules/bowtie2/align/main.nf'
include { BAMTOBEDGRAPH } from './modules/local/bedtools/main.nf'
include { SEACR_CALLPEAK } from './modules/nf-core/modules/seacr/callpeak/main.nf'
include { MACS2_CALLPEAK } from './modules/nf-core/modules/macs2/callpeak/main.nf'

//Define stdout message for the command line use
log.info """\
         C U T & R U N-  P I P E L I N E
         ===================================
         Project           : $workflow.projectDir
         Project workDir   : $workflow.workDir
         Container Engine  : $workflow.containerEngine
         Samples           : ${params.sample_sheet}
         Genome            : ${params.fasta}
         """
         .stripIndent()

//run the workflow for alignment, to bedgraph, to peak calling for Cut&Run data
//QUESTION: need a switch of MACS vs SEACR peak caller to be used?
workflow call_peaks {
    //Empty channel to collect the versions of software used 
    Channel.empty()
        .set { versions }
    //Create the input channel which contains the SAMPLE_ID, whether its single-end, and the file paths for the fastqs. 
    Channel.fromPath(file(params.sample_sheet))
        .ifEmpty { error  "No file found ${params.sample_sheet}." }
        .splitCsv(header: true, sep: '\t')
        .map { meta -> [ [ "id":meta["sample_id"], "single_end":meta["single_end"].toBoolean(), "group":meta["target_or_control"] ], //meta
                         [ file(meta["read1"], checkIfExists: true), file(meta["read2"], checkIfExists: true) ] //reads
                    ]}
        .set { meta_ch }
    //Stage th file for bedgraph generations
    Channel.fromPath(params.genome_file)
        .ifEmpty { error  "No file found ${params.genome_file}." }
        .collect() //collect converts this to a value channel and used multiple times
        .set { genome_file }
    //Stage the genome index directory
    Channel.fromPath(params.index)
        .ifEmpty { error "No directory found ${params.index}." }
        .collect() //collect converts this to a value channel and used multiple times
        .set { index }
    //Adapter and Quality trimming of the fastq files 
    TRIMGALORE(meta_ch)
    //Perform the alignement 
    BOWTIE2_ALIGN(TRIMGALORE.out.reads, index,
                  params.save_unaligned, params.sort_bam)
    //Conver the bam files to bedtools 
    // Should there be a subworkflow THE SEACR+BEDGRAPH GENERATION step? 
    BAMTOBEDGRAPH(BOWTIE2_ALIGN.out.bam, genome_file)
    BAMTOBEDGRAPH.out.bedgraph
        .view()
    /*
    //SEACR peak calling
    //Need to have channel with a tuple [meta, signal bdg, control bdg]
    BAMTOBEDGRAPH.out.bedgraph
        .map { meta -> [meta, [signal_bedgraph], [control_bedgraph] ] }
        .set { begraphs_ch }
    //SEACR_CALLPEAK(begraphs_ch)
    //MACS2 peak calling 
    BOWTIE2_ALIGN.out.bam
        .map { [meta, [ipbam], [controlbam] ] }
        .set { bam_ch }
    MACS2_CALLPEAK(bam_ch)
    */
    versions.concat(TRIMGALORE.out.versions, 
                    BOWTIE2_ALIGN.out.versions,
                    BAMTOBEDGRAPH.out.versions)
            .subscribe onNext: { println "Version: " + it }
}

//Generate the index file 
workflow bowtie2_index {
    //Stage the fasta files
    Channel.fromPath(params.fasta)
        .ifEmpty { error "No files found ${params.fasta}." }
        .set { fasta }
    //execute the BOWTIE2 genome index process
    BOWTIE2_BUILD(fasta)
    //execute the full peak calling workflow
    call_peaks()
}

//End with a message to print to standard out on workflow completion. 
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

