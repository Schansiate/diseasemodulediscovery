//
// Applies context-specific (e.g. GTEx tissue) and CRAPome filtering to graph-tool
// parsed networks, and folds the resulting filtered networks/seeds back into the
// main network/seeds channels.
//

include { CONTEXT_SPECIFIC_FILTERING                        } from '../../../modules/local/context_specific_filtering/main'
include { CONTEXT_SPECIFIC_FILTERING as CRAPOME_FILTERING   } from '../../../modules/local/context_specific_filtering/main'

workflow NETWORK_FILTERING {

    take:
    ch_network_gt 
    ch_context_specific_input  // channel: [ path(seeds), path(network), val(context), val(source), val(threshold) ]
    ch_seeds
    ch_network_multiqc
    ch_perturbed_networks

    main:
    ch_versions = channel.empty()
    ch_filtered_networks = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_node_degrees = channel.empty()
    ch_filtered_network_multiqc = channel.empty()

    ch_context_specific_network = ch_context_specific_input
        .map{_seeds, network, context_key, context, source, threshold ->
            [network.baseName, context_key, context, source, threshold]
        }
        .combine(
            ch_network_gt.map{meta, network -> [meta.network_id, meta, network]},
            by: 0
        )
        .map{ _network_id, context_key, context, source, threshold, meta, network ->
            def dup = meta.clone()
            def filtering_file = downloadFilteringFile(context, source)
            dup.id = meta.id + "." + context_key
            dup.network_id = meta.network_id + "." + context_key
            [dup, network, context, source, filtering_file, threshold]
        }
    ch_context_specific_seeds = ch_context_specific_input
        .map{seeds, network, context_key, context, source, threshold ->
            def context_specific_id = seeds.baseName + network.baseName + "." + context_key
            def network_id = network.baseName + "." + context_key
            [[id: context_specific_id, seeds_id: seeds.baseName, network_id: network_id ], seeds]
        }

    CONTEXT_SPECIFIC_FILTERING(ch_context_specific_network)
    ch_versions = ch_versions.mix(CONTEXT_SPECIFIC_FILTERING.out.versions)
    ch_filtered_networks = CONTEXT_SPECIFIC_FILTERING.out.filtered_network
    ch_filtered_network_multiqc = CONTEXT_SPECIFIC_FILTERING.out.multiqc
    ch_multiqc_files = CONTEXT_SPECIFIC_FILTERING.out.filtering_statistic
        .map{ _meta, path -> path }
        .collectFile(
            cache: false,
            storeDir: "${params.outdir}/mqc_summaries",
            name: 'filtering_statistics_mqc.tsv',
            keepHeader: true
        )
    ch_node_degrees = CONTEXT_SPECIFIC_FILTERING.out.node_degree
    
    if(params.filter_crapomes){
        CRAPOME_FILTERING(
            ch_network_gt.mix(ch_filtered_networks)
            .map{meta, network ->
                def dup = meta.clone()
                def filtering_file = downloadFilteringFile("","CRAPOME")
                dup.id = meta.id + ".crapome" + "." + params.crapome_filtering_threshold
                dup.network_id = meta.network_id + ".crapome"+"."+ params.crapome_filtering_threshold
                [dup, network, "CRAPome", "CRAPome", filtering_file, params.crapome_filtering_threshold]
            }
        )
        ch_filtered_networks = ch_filtered_networks.mix(CRAPOME_FILTERING.out.filtered_network)
        ch_filtered_network_multiqc = ch_filtered_network_multiqc.mix(CRAPOME_FILTERING.out.multiqc)
        ch_multiqc_files = ch_multiqc_files
        .mix(CRAPOME_FILTERING.out.filtering_statistic
            .map{ _meta, path -> path }
            .collectFile(
                cache: false,
                storeDir: "${params.outdir}/mqc_summaries",
                name: 'crapome_filtering_statistics_mqc.tsv',
                keepHeader: true
            ))
         ch_crapome_specific_seeds = ch_seeds.mix(ch_context_specific_seeds)
            .map{meta, seeds ->
                def dup = meta.clone()
                dup.id = meta.id + ".crapome" + "." + params.crapome_filtering_threshold
                dup.network_id = meta.network_id + ".crapome" + "." + params.crapome_filtering_threshold
                [dup, seeds]
            }
        ch_context_specific_seeds = ch_context_specific_seeds.mix(ch_crapome_specific_seeds)
        ch_node_degrees = ch_node_degrees.mix(CRAPOME_FILTERING.out.node_degree)

    }
    if(params.run_filtered_networks_only){
        ch_network_gt = ch_filtered_networks
        ch_seeds = ch_context_specific_seeds
        ch_network_multiqc = ch_filtered_network_multiqc
        ch_perturbed_networks = ch_filtered_networks.map{ meta, _path -> [meta, []] }
    } else {
        ch_network_gt = ch_network_gt.mix(ch_filtered_networks)
        ch_seeds = ch_seeds.mix(ch_context_specific_seeds)
        ch_network_multiqc = ch_network_multiqc.mix(ch_filtered_network_multiqc)
        ch_perturbed_networks = ch_perturbed_networks.mix(ch_filtered_networks.map{ meta, _path -> [meta, []] })
    }

    emit:
    versions              = ch_versions                            // channel: [ path(versions.yml) ]
    network_gt             = ch_network_gt                          // channel: [ val(meta[id,network_id]), path(network) ]
    seeds                  = ch_seeds                                // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    network_multiqc        = ch_network_multiqc                     // channel: path(multiqc)
    perturbed_networks     = ch_perturbed_networks                  // channel: [ val(meta[id,network_id]), [path(perturbed_networks)] ]
    multiqc_files          = ch_multiqc_files
    expression_distribution = CONTEXT_SPECIFIC_FILTERING.out.expression_distribution
    node_degree             = ch_node_degrees        // channel: [ val(meta), path(node_degree) ]
   
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FILTERING FILE DOWNLOAD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def downloadFilteringFile(context, source){
    if (source == "GTEX"){
        return file("https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz")
    }
    else if (source == "PAXDB"){
        return file("https://pax-db.org/downloads/6.1/datasets/9606/9606-${context}-integrated.txt")
    }
    else if (source == "TCGA"){
        return file("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-${context}.star_tpm.tsv.gz")
    }
    else if (source == "CRAPOME"){
        return file("https://reprint-apms.org/?q=system/files/crap_db_v1_flat_file_human.xlsx")
    }
}