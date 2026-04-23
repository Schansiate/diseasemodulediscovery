#!/usr/bin/env python

import argparse
from pathlib import Path
import sys
import pandas as pd 
import graph_tool.all as gt 
import util as utils 

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="filter network by tissue specific expression")
    parser.add_argument("--network",
                        type=str,
                        required=True,
                        help="path to the network file in gt format")
    parser.add_argument("--tissues",
                        type=lambda x: x.split(','),
                        required=True,
                        help="tissues to filter by")
    parser.add_argument("--threshold",
                        type=float,
                        default=1.0,
                        help="expression threshold to filter by")
    parser.add_argument("--expression_file",
                        type=str,
                        required=True,
                        help="path to the expression file")
    return parser.parse_args(argv)

def filter_network_by_tissue(network_file, tissues, expression_file, threshold):
    network = gt.load_graph(network_file)
    stem = Path(network_file).stem
    name_index = utils.name2index(network)
    expression_by_tissue = pd.read_csv(expression_file, sep="\t", skiprows=2, header = 0)
    expression_by_tissue["Name"] = expression_by_tissue["Name"].apply(lambda x: x.split(".")[0])
    in_network_genes = expression_by_tissue[expression_by_tissue["Name"].isin(name_index.keys())].copy()
    in_network_genes["vertex_id"] = in_network_genes["Name"].map(name_index)
    for tissue in tissues:
        tissue_specific_genes = in_network_genes[in_network_genes[tissue] > threshold]
        tissue_filter = network.new_vertex_property("bool")
        tissue_filter.a[tissue_specific_genes["vertex_id"]] = True
        network.set_vertex_filter(tissue_filter)
        tissue_specific_graph = gt.GraphView(network, vfilt=tissue_filter)
        tissue_specific_graph.save(f"{stem}_{tissue}_specific.gt")



def main(argv=None):
    args = parse_args(argv)
    filter_network_by_tissue(args.network ,args.tissues, args.expression_file, args.threshold)

if __name__ == "__main__":
    sys.exit(main())
