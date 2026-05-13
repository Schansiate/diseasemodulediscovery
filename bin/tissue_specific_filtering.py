#!/usr/bin/env python

import argparse
from pathlib import Path
import sys
import pandas as pd 
import graph_tool.all as gt 
import util as utils 
from gprofiler import GProfiler 

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="filter network by tissue specific expression")
    parser.add_argument("--network",
                        type=str,
                        required=True,
                        help="path to the network file in gt format")
    parser.add_argument("--tissue",
                        required=True,
                        help="tissue to filter by")
    parser.add_argument("--threshold",
                        type=float,
                        default=1.0,
                        help="expression threshold to filter by")
    parser.add_argument("--expression_file",
                        type=str,
                        required=True,
                        help="path to the expression file")
    parser.add_argument("--id_space",
                        type=str,
                        default="ensembl",
                        help="id space of the gene names in the network and expression file")
    return parser.parse_args(argv)


def convert_id_space(expression_df, id_space):
    ids = expression_df["Name"]
    if id_space == "ensembl": 
        return expression_df
    elif id_space == "entrez":
        targetSpace = "ENTREZGENE_ACC"
    elif id_space == "uniprot":
        targetSpace = "UNIPROT_GN_ACC"
    else:
        raise ValueError(f"Unsupported id_space: {id_space}")
    gp = GProfiler(return_dataframe=True)
    query_result = gp.convert(organism="hsapiens", query=ids.tolist(),target_namespace=targetSpace)
    #not_found = query_result.loc[query_result["converted"].astype(str) == "None", "incoming"].shape[0] / expression_df.shape[0]
    query_result = query_result.merge(expression_df, left_on="incoming", right_on="Name")
    query_result = query_result[["incoming", "converted", "expression"]]
    collapsed_result = (
        query_result[query_result["converted"].astype(str) != "None"]
        .assign(converted=lambda x : x["converted"].astype(str))
        .drop(columns=["incoming"])
        .groupby("converted", as_index=False)
        .agg({
            "expression": "sum"
        })
    )
    collapsed_result = collapsed_result.rename(columns={"converted": "Name"})
    return collapsed_result



def filter_network(network_file, threshold, expression_by_tissue, tissue):
    network = gt.load_graph(network_file)
    stem = Path(network_file).stem
    name_index = utils.name2index(network)
    tissue_specific_genes = expression_by_tissue[expression_by_tissue["expression"] > threshold]
    in_network_tissue_genes = tissue_specific_genes[tissue_specific_genes["Name"].isin(name_index.keys())].copy()
    not_in_network = tissue_specific_genes.loc[~tissue_specific_genes["Name"].isin(name_index.keys())].shape[0]
    in_network_tissue_genes["vertex_id"] = in_network_tissue_genes["Name"].map(name_index)
    tissue_filter = f"{tissue}_filter"
    network.vp[tissue_filter] = network.new_vertex_property("bool")
    for vertex_id in in_network_tissue_genes["vertex_id"].unique():
        network.vp[tissue_filter][vertex_id] = True
    in_Network_filtered_genes = network.num_vertices() - len(in_network_tissue_genes["vertex_id"].unique())
    network.set_vertex_filter(network.vp[tissue_filter])
    network.purge_vertices()
    network.clear_filters()
    del network.vp[tissue_filter]
    network.save(f"{stem}.{tissue}.gt")
    
    return tissue_specific_genes.shape[0], not_in_network, in_Network_filtered_genes



def main(argv=None):
    args = parse_args(argv)
    global id_space 
    id_space = args.id_space
    expression_by_tissue = pd.read_csv(args.expression_file, sep="\t", skiprows=2, header=0)
    #trim version numbers from ensembl IDs 
    expression_by_tissue["Name"] = expression_by_tissue["Name"].apply(lambda x: x.split(".")[0])
    #select only the relevant tissue and name columns
    expression_by_tissue = expression_by_tissue[["Name", args.tissue]]
    expression_by_tissue = expression_by_tissue.rename(columns={args.tissue: "expression"})
    #map gene_ids to the specified id space 
    expression_by_tissue = convert_id_space(expression_by_tissue, args.id_space)
    #filter network by tissue specific expression with the given threshold 
    n_genes_in_tissue, n_not_in_network, in_Network_filtered_genes = filter_network(args.network ,args.threshold, expression_by_tissue, args.tissue)
    with open("filtering_statistic.tsv", "w") as f:
        f.write("tissue\tgenes_in_tissue\tgenes_not_in_network\tin_network_filtered_genes\n")
        f.write(f"<a href=https://gtexportal.org/home/tissue/{args.tissue}>{args.tissue}</a>\t{n_genes_in_tissue}\t{n_not_in_network}\t{in_Network_filtered_genes}\n")
if __name__ == "__main__":
    sys.exit(main())
