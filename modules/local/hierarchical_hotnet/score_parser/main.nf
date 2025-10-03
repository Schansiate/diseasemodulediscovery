process HIERARCHICAL_HOTNET_SCORE_PARSER {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta),  path(node_list)
    output:
    tuple val(meta), path("*.node_scores.tsv")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat "/Users/motan/Research/tum/diseasemodulediscovery/work/0e/4c21176797e8316a49bd88deea8cef/scores_1.tsv" > "${meta.id}.node_scores.tsv"
    """
}
