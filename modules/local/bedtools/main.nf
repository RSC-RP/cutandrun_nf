// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.
process BAMTOBEDGRAPH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3':
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(bam)
    path chrom_sizes
    val spike_norm
    val scale_factor

    output:
    tuple val(meta), path("*fragments.bg"),    emit: bedgraph
    tuple val(meta), path("*.bed"),             emit: bedfiles
    path "versions.yml"           ,             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def scale_arg = spike_norm ? "-scale $scale_factor" : ''
    def suffix = spike_norm ? '_normalized' : ''
    """
    set -eu -o pipefail

    bedtools bamtobed \\
        -bedpe \\
        -i ${bam} \\
        $args > ${prefix}_aligned.bed

    # Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' ${prefix}_aligned.bed > ${prefix}_aligned.clean.bed
    cut -f 1,2,6 ${prefix}_aligned.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${prefix}_aligned.fragments.bed
    
    if [[ $spike_norm ]]
    then
        echo "The spikein scale factor is: $scale_factor"
    else

    bedtools genomecov \\
        $scale_arg \\
        $args \\
        -bg \\
        -i ${prefix}_aligned.fragments.bed \\
        -g ${chrom_sizes} > ${prefix}_aligned${suffix}_fragments.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
