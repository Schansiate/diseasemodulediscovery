process TISSUE_SPECIFIC_FILTERING {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(network)
    path(expression_file)
    val tissues

    output:
    tuple val(meta), path("${meta.id}_*specific.gt"), emit: filtered_networks
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    tissue_specific_filtering.py --network ${network} --tissues ${tissues} --expression_file ${expression_file}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}