nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from './modules/nf-core/modules/bowtie2/build/main.nf'
include { TRIMGALORE } from './modules/nf-core/modules/trimgalore/main.nf'
include { BOWTIE2_ALIGN } from './modules/nf-core/modules/bowtie2/align/main.nf'
include { BAMTOBEDGRAPH } from './modules/local/bedtools/main.nf'
include { SEACR_CALLPEAK } from './modules/nf-core/modules/seacr/callpeak/main.nf'
include { MACS2_CALLPEAK } from './modules/nf-core/modules/macs2/callpeak/main.nf'
include { KHMER_UNIQUEKMERS } from './modules/nf-core/modules/khmer/uniquekmers/main.nf'

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


/*
https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
https://khmer.readthedocs.io/en/stable/user/scripts.html#unique-kmers-py
*/

//run the workflow for alignment, to bedgraph, to peak calling for Cut&Run data
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
        // Conver the bam files to bed format with bedtools 
        BAMTOBEDGRAPH(BOWTIE2_ALIGN.out.bam, genome_file)
        //Define the control and the target channels. should be a custom groovy function really - need to figure this out.
        BAMTOBEDGRAPH.out.bedgraph
            .branch { 
               control: it[0].group =~ /control/
               targets: it[0].group =~ /target/
            }
            .set { bedgraphs }
        bedgraphs.control
            .cross(bedgraphs.targets){ meta -> meta[0].sample } // join by the key name "sample"
            .map { meta -> [ meta[1][0], meta[1][1], meta[0][1] ] }
            .set { seacr_ch }
        //SEACR peak calling
        SEACR_CALLPEAK(seacr_ch, params.threshold)

        // MACS2 peak calling , Optionally run macs2 peak calling
        // can the channel creation be done in the subworkflow?
        if ( params.run_macs2 ){
            BOWTIE2_ALIGN.out.bam
                .branch { 
                    control: it[0].group =~ /control/
                    targets: it[0].group =~ /target/
                }
                .set { bams }
            bams.control
                .cross(bams.targets){ meta -> meta[0].sample }
                .map { meta -> [ meta[1][0], meta[1][1], meta[0][1] ] }
                .set { macs_ch }
            //Run MAC2 peak calling
            macs2_peaks(macs_ch)
        }
        // versions.concat(TRIMGALORE.out.versions, 
        //                 BOWTIE2_ALIGN.out.versions,
        //                 BAMTOBEDGRAPH.out.versions, 
        //                 SEACR_CALLPEAK.out.versions)
        //         .collect()
        //         .collectFile(name: 'versions.txt', newLine: true)
}

//Generate the index file (subworkflow)
workflow bowtie2_index {
    main:        
    //Stage the fasta files
    Channel.fromPath(file(params.fasta, checkIfExists: true))
        .set { fasta }
    //execute the BOWTIE2 genome index process
    BOWTIE2_BUILD(fasta)

    emit:
    index = BOWTIE2_BUILD.out.index
    versions = BOWTIE2_BUILD.out.versions
}

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


//End with a message to print to standard out on workflow completion. 
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

