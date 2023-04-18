// Function to determine which container should be used in the workflow based on SEACR version requested
def select_container (VERSION){
    containers = [ sif:'', docker:'' ] // assoc array
    if ( VERSION == '1.3' ){
        sif = 'https://depot.galaxyproject.org/singularity/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:f4bb19b68e66de27e4c64306f951d5ff11919931-0'
        docker = 'quay.io/biocontainers/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:f4bb19b68e66de27e4c64306f951d5ff11919931-0'
    }else if ( VERSION == '1.4' ) {
        sif = 'docker://quay.io/jennylsmith/seacr:v1.4'
        docker = 'quay.io/jennylsmith/seacr:v1.4'
    }
    containers.sif = sif 
    containers.docker = docker
    return(containers)
}

process SEACR_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::seacr=1.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? select_container(task.ext.version).sif : select_container(task.ext.version).docker }"

    input:
    tuple val(meta), path(bedgraph), path(ctrlbedgraph)
    val threshold
    val spike_norm

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = "${task.ext.version}" // Version information not provided by tool on CLI
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def function_switch = ctrlbedgraph ? "$ctrlbedgraph" : "$threshold"
    def spikein = spike_norm ?  'spikein_norm' : 'non'
    def type = "${function_switch.toString().replaceAll("_aligned.+","")}"
    def method =  "${args}" =~ /\bnorm\b/ ? 'norm' : "${spikein}"
    def suffix = ctrlbedgraph ? "vs_${type}_${method}" : "threshold${type}_${method}"
    if ( VERSION == '1.4' )
        """
        SEACR_1.4.sh \\
            -b $bedgraph \\
            -c $function_switch \\
            $args \\
            -o ${prefix}_${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seacr: $VERSION
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        END_VERSIONS
        """
    else if ( VERSION == '1.3' )
        """
        SEACR_1.3.sh \\
            $bedgraph \\
            $function_switch \\
            $args \\
            ${prefix}_${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seacr: $VERSION
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        END_VERSIONS
        """
    else
        error "Invalid SEACR version: ${VERSION}. Please choose 1.3 or 1.4"
}
