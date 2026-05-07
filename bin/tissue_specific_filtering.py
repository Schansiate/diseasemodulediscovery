#!/usr/bin/env python

import argparse
from pathlib import Path
import sys
import pandas as pd 
import graph_tool.all as gt 
import util as utils 
from gprofiler import GProfiler 
import logging 

logger = logging.getLogger()

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
    #TODO: add loggin of not found ids 
    not_found = query_result.loc[query_result["converted"].astype(str) == "None", "incoming"].nunique()
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
    in_network_tissue_genes["vertex_id"] = in_network_tissue_genes["Name"].map(name_index)
    tissue_filter = network.new_vertex_property("bool")
    tissue_filter.a[in_network_tissue_genes["vertex_id"]] = True
    network.set_vertex_filter(tissue_filter)
    tissue_specific_graph = gt.GraphView(network, vfilt=tissue_filter)
    tissue_specific_graph.save(f"{stem}.{tissue}.gt")



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
    print(expression_by_tissue.head())
    filter_network(args.network ,args.threshold, expression_by_tissue, args.tissue)

if __name__ == "__main__":
    sys.exit(main())
