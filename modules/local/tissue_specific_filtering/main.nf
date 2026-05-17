process TISSUE_SPECIFIC_FILTERING {
    tag "$meta.id"
    label 'process_single'
    container 'docker.io/motan04/modulediscovery_python_dependencies:latest'
    input:
    tuple val(meta), path(network), val(tissue)
    path(expression_file)

    output:
    tuple val(meta), path("${meta.id}.gt")             , emit: filtered_network
    tuple val(meta), path("input_network_multiqc.tsv") , emit: multiqc
    tuple val(meta), path("filtering_statistic.tsv")  , emit: filtering_statistic
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tissue_specific_filtering.py --network ${network} \
    --tissue ${tissue} \
    --expression_file ${expression_file} \
    --id_space ${params.id_space} \
    --threshold ${params.filtering_threshold}

    graph_tool_parser.py ${meta.id}.gt -f gt -l DEBUG
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
