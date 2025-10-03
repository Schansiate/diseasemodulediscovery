process HIERARCHICAL_HOTNET_INPUT_PARSER {
    tag "meta.id"
    label 'process_single'
    
    input:
    tuple val(meta), (path(network))
    
    output:
    tuple val(meta), path("*.node_list.tsv"), path("*.edge_list.tsv") , emit: network
    path "versions.yml"                                         , emit: versions 

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    cat "/Users/motan/Research/tum/diseasemodulediscovery/work/0e/4c21176797e8316a49bd88deea8cef/network_1_edge_list.tsv" > "${meta.id}.edge_list.tsv"
    cat "/Users/motan/Research/tum/diseasemodulediscovery/work/0e/4c21176797e8316a49bd88deea8cef/network_1_index_gene.tsv" > "${meta.id}.node_list.tsv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}