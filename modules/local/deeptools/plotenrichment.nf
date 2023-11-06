process DEEPTOOLS_PLOTENRICHMENT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::deeptools=3.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)

    output:
    tuple val(meta), path("*.txt"), emit: stats
    tuple val(meta), path("*.pdf"), emit: plot
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotEnrichment \\
        $args \\
        -b $bam \\
        --BED $bed \\
        -o ${prefix}.plotEnrichment.pdf \\
        --outRawCounts ${prefix}.plotEnrichment.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotEnrichment --version | sed -e "s/plotEnrichment //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf
    touch ${prefix}.txt 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotPCA --version | sed -e "s/plotPCA //g")
    END_VERSIONS
    """
}
