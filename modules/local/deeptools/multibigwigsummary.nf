process DEEPTOOLS_MULTIBIGWIGSUMMARY {
    tag "$meta.id"
    label 'process_medium'

     conda "bioconda::deeptools=3.5.2"
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0':
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0' }"


    input:
    tuple val(meta), path(bigwigs)

    output:
    tuple val(meta), path("*.npz"), emit: npz
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    multiBigwigSummary bins \
        -b ${bigwigs} \
        -out ${prefix}_scores_per_bin.npz \
        --outRawCounts ${prefix}_scores_per_bin.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBigwigSummary --version | sed -e "s/multiBigwigSummary //g")
    END_VERSIONS
    """
}
