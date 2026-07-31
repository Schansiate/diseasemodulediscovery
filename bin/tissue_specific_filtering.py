#!/usr/bin/env python

import argparse
from pathlib import Path
import sys
import requests
import pandas as pd
import yaml
import graph_tool.all as gt
import util as utils
from gprofiler import GProfiler

GTEX_URL = "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
PAXDB_URL = (
    "https://pax-db.org/downloads/6.1/datasets/9606/9606-{tissue}-integrated.txt"
)
TCGA_URL = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{cancer_type}.star_tpm.tsv.gz"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="filter network by tissue specific expression"
    )
    parser.add_argument(
        "--network",
        type=str,
        required=True,
        help="path to the network file in gt format",
    )
    parser.add_argument("--tissue", required=True, help="tissue to filter by")
    parser.add_argument(
        "--threshold", type=float, default=1.0, help="expression threshold to filter by"
    )
    parser.add_argument(
        "--filtering_source",
        type=str,
        default=None,
        help="database name to use as expression source (GTEx, ProteomicsDB, PAXDB, TCGA)",
    )
    parser.add_argument(
        "--custom_expression_file",
        type=str,
        default=None,
        help="path to a custom expression file with columns 'id' and 'expression'",
    )
    parser.add_argument(
        "--id_space",
        type=str,
        default="ensembl",
        help="id space of the gene names in the network and expression file",
    )
    return parser.parse_args(argv)


def resolve_expression_file(filtering_source, custom_expression_file, tissue):
    """Returns a tuple of (expression_df, dataset_link). dataset_link is a URL
    pointing at the specific dataset used for filtering, or None if the
    source has no such page (custom files)."""
    if filtering_source != "null" and custom_expression_file != "null":
        raise ValueError(
            "Specify either --filtering_source or --custom_expression_file, not both"
        )
    if filtering_source == "null" and custom_expression_file == "null":
        raise ValueError(
            "One of --filtering_source or --custom_expression_file must be specified"
        )

    if custom_expression_file != "null":
        df = pd.read_csv(custom_expression_file, sep="\t", header=0)
        df = df[["id", tissue]].rename(columns={tissue: "expression"})
        return df, None

    if filtering_source == "GTEx":
        response = requests.get(GTEX_URL, stream=True)
        response.raise_for_status()
        # write to the task work directory (rather than a python tempfile) so
        # the raw source data is kept alongside the rest of the process outputs
        expression_file = Path("GTEx_expression.gct.gz")
        with open(expression_file, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        df = pd.read_csv(expression_file, sep="\t", skiprows=2, header=0)
        df.rename(columns={df.columns[0]: "id"}, inplace=True)
        df = df[["id", tissue]].rename(columns={tissue: "expression"})
        df["id"] = df["id"].apply(lambda x: x.split(".")[0])
        return df, f"https://gtexportal.org/home/tissue/{tissue}"
    elif filtering_source == "PAXDB":
        response = requests.get(PAXDB_URL.format(tissue=tissue.upper()))
        response.raise_for_status()
        expression_file = Path(f"9606-{tissue}-integrated.txt")
        expression_file.write_text(response.text)
        # the first line of the file holds the PaxDB dataset id, e.g. "#id: 4032974247"
        dataset_id = response.text.splitlines()[0].split(":", 1)[1].strip()
        # the column header itself is also '#'-prefixed in PaxDB files, so it
        # gets dropped as a comment too; supply the column names explicitly
        df = pd.read_csv(
            expression_file,
            sep="\t",
            comment="#",
            header=None,
            names=["gene_name", "string_external_id", "abundance"],
        )
        df = df.rename(columns={"gene_name": "id", "abundance": "expression"})
        df = df[["id", "expression"]]
        return df, f"https://pax-db.org/dataset/9606/{dataset_id}/"
    elif filtering_source == "TCGA":
        if not tissue.upper().startswith("TCGA_"):
            raise ValueError(
                f"--tissue must be of the form 'TCGA_{{cancer_type}}' when using --filtering_source TCGA, got '{tissue}'"
            )
        cancer_type = tissue.upper().removeprefix("TCGA_")
        response = requests.get(TCGA_URL.format(cancer_type=cancer_type), stream=True)
        response.raise_for_status()
        expression_file = Path(f"TCGA-{cancer_type}.star_tpm.tsv.gz")
        with open(expression_file, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        df = pd.read_csv(expression_file, sep="\t", header=0)
        df.rename(columns={df.columns[0]: "id"}, inplace=True)
        df["id"] = df["id"].apply(lambda x: x.split(".")[0])
        sample_columns = df.columns.drop("id")
        # values are log2(tpm+1) transformed; undo that before taking the median
        # so it is computed on the same TPM scale as GTEx
        tpm = 2 ** df[sample_columns] - 1
        df["expression"] = tpm.median(axis=1)
        df = df[["id", "expression"]]
        return df, None
    elif filtering_source == "ProteomicsDB":
        raise NotImplementedError(f"{filtering_source} support is not yet implemented")
    else:
        raise ValueError(
            f"Unknown filtering_source: '{filtering_source}'. Must be one of: GTEx, ProteomicsDB, PAXDB, TCGA"
        )


def convert_id_space(expression_df, id_space):
    ids = expression_df["id"]
    if id_space == "ensembl":
        return expression_df
    elif id_space == "entrez":
        targetSpace = "ENTREZGENE_ACC"
    elif id_space == "uniprot":
        targetSpace = "UNIPROT_GN_ACC"
    elif id_space == "symbol":
        targetSpace = "HGNC"
    else:
        raise ValueError(f"Unsupported id_space: {id_space}")
    gp = GProfiler(return_dataframe=True)
    query_result = gp.convert(
        organism="hsapiens", query=ids.tolist(), target_namespace=targetSpace
    )
    # not_found = query_result.loc[query_result["converted"].astype(str) == "None", "incoming"].shape[0] / expression_df.shape[0]
    query_result = query_result.merge(expression_df, left_on="incoming", right_on="id")
    query_result = query_result[["incoming", "converted", "expression"]]
    collapsed_result = (
        query_result[query_result["converted"].astype(str) != "None"]
        .assign(converted=lambda x: x["converted"].astype(str))
        .drop(columns=["incoming"])
        .groupby("converted", as_index=False)
        .agg({"expression": "sum"})
    )
    collapsed_result = collapsed_result.rename(columns={"converted": "id"})
    return collapsed_result


def save_expression_distribution(expression_by_tissue, stem, tissue, threshold):
    values = expression_by_tissue["expression"]
    if values.empty:
        return
    percentiles = list(range(0, 101))
    pct_values = [float(values.quantile(p / 100)) for p in percentiles]
    distribution = {
        "name": f"{stem}.{tissue}",
        "threshold": threshold,
        "data": [[v, p] for p, v in zip(percentiles, pct_values)],
    }
    with open(f"{stem}.{tissue}.expression_distribution.yaml", "w") as f:
        yaml.safe_dump(distribution, f, sort_keys=False, default_flow_style=None)


def filter_network(network_file, threshold, expression_by_tissue, tissue):
    network = gt.load_graph(network_file)
    stem = Path(network_file).stem
    name_index = utils.name2index(network)
    network_num_vertices = network.num_vertices()
    expression_num_entries = expression_by_tissue.shape[0]
    in_network_expression = expression_by_tissue[
        expression_by_tissue["id"].isin(name_index.keys())
    ].copy()
    genes_not_in_network = expression_num_entries - in_network_expression.shape[0]
    genes_not_in_expression_file = (
        network_num_vertices - in_network_expression["id"].nunique()
    )
    tissue_specific_in_network = in_network_expression[
        in_network_expression["expression"] > threshold
    ].copy()
    genes_filtered_by_threshold = (
        in_network_expression.shape[0] - tissue_specific_in_network.shape[0]
    )
    tissue_specific_in_network["vertex_id"] = tissue_specific_in_network["id"].map(
        name_index
    )

    # filter network by tissue specific expression
    network.vp["tissue_filter"] = network.new_vertex_property("bool")
    network.vp["expression_in_tissue"] = network.new_vertex_property("double")
    for vertex_id in tissue_specific_in_network["vertex_id"].unique():
        network.vp["tissue_filter"][vertex_id] = True
        expression_value = tissue_specific_in_network.loc[
            tissue_specific_in_network["vertex_id"] == vertex_id
        ]["expression"].values[0]
        network.vp["expression_in_tissue"][vertex_id] = expression_value
    network.set_vertex_filter(network.vp["tissue_filter"])
    network.purge_vertices()
    network.clear_filters()
    del network.vp["tissue_filter"]
    network.save(f"{stem}.{tissue}.gt")
    save_expression_distribution(in_network_expression, stem, tissue, threshold)

    return {
        "genes_not_in_network_absolute": genes_not_in_network,
        "genes_not_in_network_relative": (
            genes_not_in_network / expression_num_entries
            if expression_num_entries
            else 0.0
        ),
        "genes_not_in_expression_file_absolute": genes_not_in_expression_file,
        "genes_not_in_expression_file_relative": (
            genes_not_in_expression_file / network_num_vertices
            if network_num_vertices
            else 0.0
        ),
        "genes_filtered_by_threshold_absolute": genes_filtered_by_threshold,
        "genes_filtered_by_threshold_relative": (
            genes_filtered_by_threshold / network_num_vertices
            if network_num_vertices
            else 0.0
        ),
    }


def main(argv=None):
    args = parse_args(argv)
    global id_space
    id_space = args.id_space
    expression_by_tissue, dataset_link = resolve_expression_file(
        args.filtering_source, args.custom_expression_file, args.tissue
    )
    # map gene_ids to the specified id space
    expression_by_tissue = convert_id_space(expression_by_tissue, args.id_space)
    n_genes_in_tissue = expression_by_tissue[
        expression_by_tissue["expression"] > args.threshold
    ].shape[0]
    # filter network by tissue specific expression with the given threshold
    filtering_statistics = filter_network(
        args.network, args.threshold, expression_by_tissue, args.tissue
    )
    tissue_label = (
        f"<a href={dataset_link}>{args.tissue}</a>" if dataset_link else args.tissue
    )
    with open("filtering_statistic.tsv", "w") as f:
        f.write(
            "network\ttissue\tgenes_in_tissue\t"
            "genes_not_in_network_absolute\tgenes_not_in_network_relative\t"
            "genes_not_in_expression_file_absolute\tgenes_not_in_expression_file_relative\t"
            "genes_filtered_by_threshold_absolute\tgenes_filtered_by_threshold_relative\n"
        )
        f.write(
            f"{Path(args.network).stem}.{args.tissue}\t{tissue_label}\t{n_genes_in_tissue}\t"
            f"{filtering_statistics['genes_not_in_network_absolute']}\t{filtering_statistics['genes_not_in_network_relative']}\t"
            f"{filtering_statistics['genes_not_in_expression_file_absolute']}\t{filtering_statistics['genes_not_in_expression_file_relative']}\t"
            f"{filtering_statistics['genes_filtered_by_threshold_absolute']}\t{filtering_statistics['genes_filtered_by_threshold_relative']}\n"
        )


if __name__ == "__main__":
    sys.exit(main())
