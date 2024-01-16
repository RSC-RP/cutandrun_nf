process MACSPEAKSTOBED {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3':
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(peaks)

    output:
    tuple val(meta), path("*.bed")      , emit: bed
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f 1,2,3 $peaks > ${prefix}_peaks.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_peaks.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
