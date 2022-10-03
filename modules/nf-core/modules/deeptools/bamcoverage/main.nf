process DEEPTOOLS_BAMCOVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::deeptools=3.5.1 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:2c687053c0252667cca265c9f4118f2c205a604c-0':
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:2c687053c0252667cca265c9f4118f2c205a604c-0' }"

    input:
    tuple val(meta), path(input), path(input_index)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.bigWig")        ,  emit: bigwig, optional: true
    tuple val(meta), path("*.bedgraph")      ,  emit: bedgraph, optional: true
    path "versions.yml"                      ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // cram_input is currently not working with deeptools
    // therefore it's required to convert cram to bam first
    def is_cram = input.Extension == "cram" ? true : false
    def input_out = is_cram ? input.BaseName + ".bam" : "${input}"
    def sorted_bam = input.simpleName + ".coord.sort.bam"
    def fai_reference = fasta_fai ? "--fai-reference ${fasta_fai}" : ""
    def method = ( args =~ /--normalizeUsing (CPM|RPKM|BPM|RPGC)/ )
    def suffix = method.find() ? method.group(1) : ''
    def file_ext = args.contains('--outFileFormat bedgraph') ? 'bedgraph' : 'bigWig'
    if (is_cram){
    """
        if [[ $is_cram ]]
        then
            samtools view -T $fasta $input $fai_reference -@ $task.cpus -o $input_out
        fi

        samtools sort -o $sorted_bam $input_out
        samtools index -b $sorted_bam -@ $task.cpus

        bamCoverage \\
            --bam $sorted_bam \\
            $args \\
            --numberOfProcessors ${task.cpus} \\
            --outFileName ${prefix}_${suffix}.${file_ext}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
        END_VERSIONS
    """
    } else {
    """
    bamCoverage \\
        --bam $input_out \\
        $args \\
        --numberOfProcessors ${task.cpus} \\
        --outFileName ${prefix}_${suffix}.${file_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """
    }
}
