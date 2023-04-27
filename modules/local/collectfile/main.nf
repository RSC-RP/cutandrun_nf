process COLLECTFILE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'ubuntu:20.04' }"

    input:
    tuple val(meta), val(outfile)
    path scale_factor

    output:
    path "*.csv",             emit: scale_factors

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls -alh
    echo id,single_end,group,sample,spike_seq_depth,scale_factor > "${outfile}"
    cat "${scale_factor}" >> "${outfile}"
    """
}
