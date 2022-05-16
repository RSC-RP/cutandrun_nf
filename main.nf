nextflow.enable.dsl = 2

include { TRIMGALORE } from './modules/nf-core/modules/nf-core/trimgalore/main.nf'
include { BOWTIE2_ALIGN } from './modules/nf-core/modules/nf-core/bowtie2/align/main.nf'
include { BOWTIE2_BUILD } from './modules/nf-core/modules/nf-core/bowtie2/build/main.nf'
include { SEACR_CALLPEAK } from './modules/nf-core/modules/nf-core/seacr/callpeak/main.nf'
include { MACS2_CALLPEAK } from './modules/nf-core/modules/nf-core/macs2/callpeak.main.nf'

//Define stdout message for the command line use
log.info """\
         XXXX-  P I P E L I N E
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
    //Create the input channel which contains the SAMPLE_ID, whether its single-end, and the file paths for the fastqs. 
    Channel.fromPath(file(params.sample_sheet))
        .ifEmpty { error  "No file found ${params.sample_sheet}." }
        .splitCsv(header: true, sep: '\t')
        .map { meta -> [ [ "id":meta["id"], "single_end":meta["single_end"].toBoolean(),
                            "group":meta["target_or_control"]  ], //meta
                         [ file(meta["r1"], checkIfExists: true), file(meta["r2"], checkIfExists: true) ] //reads
                    ]}
        .set { meta_ch }
    //Stage th file for STAR aligner
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
    BAMTOBEDGRAPH(BOWTIE2_ALIGN.out.bam, genome_file)
    //SEACR peak calling - *SHOULD BE A SUBWORKFLOW FOR THE BEDGRAPH GENERATION??*
    //Need to have channel with a tuple [meta, signal bdg, control bdg]
    BAMTOBEDGRAPH.out.bedgraph
        .map { meta -> [meta, [signal_bedgraph], [control_bedgraph] ] }
        .set { begraphs_ch }
    SEACR_CALLPEAK(begraphs_ch)

}

//Generate the index file 
workflow bowtie_index {
    //Stage the fasta files
    Channel.fromPath(params.fasta)
        .ifEmpty { error "No files found ${params.fasta}." }
        .set {fasta}
    //execute the BOWTIE2 genome index process
    BOWTIE2_BUILD(fasta)
}

//End with a message to print to standard out on workflow completion. 
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

