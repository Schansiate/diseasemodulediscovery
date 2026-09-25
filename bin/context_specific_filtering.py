#!/usr/bin/env python

import argparse
from pathlib import Path
import sys
import pandas as pd
import yaml
import graph_tool.all as gt
import util as utils
from gprofiler import GProfiler


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
     
    )
    parser.add_argument(
        "--network",
        type=str,
        required=True,
        help="path to the network file in gt format",
    )
    parser.add_argument("--context", required=True, help="context to filter by")
    parser.add_argument(
        "--threshold",
        type=str,
        default="0.1",
        help="expression threshold to filter by. For --filtering_source CRAPome, "
        "this is the fraction (0-1) of control experiments a protein may be "
        "detected in before it is considered a contaminant",
    )
    parser.add_argument(
        "--filtering_source",
        type=str,
        default=None,
        help="database name to use as expression source (GTEx, ProteomicsDB, PAXDB, TCGA, CRAPome)",
    )
    parser.add_argument(
        "--custom_expression_file",
        type=str,
        default=None,
        help="path to a custom expression file with columns 'id' and 'expression'",
    )
    parser.add_argument(
        "--filtering_file",
        type=str,
        default=None,
        help="path to the filtering source file for --filtering_source, pre-downloaded by Nextflow",
    )
    parser.add_argument(
        "--id_space",
        type=str,
        default="ensembl",
        help="id space of the gene names in the network and expression file",
    )
    return parser.parse_args(argv)


def resolve_expression_file(
    filtering_source, custom_expression_file, filtering_file, context
):
    """Returns a tuple of (expression_df, dataset_link). dataset_link is a URL
    pointing at the specific dataset used for filtering, or None if the
    source has no such page (custom files). filtering_file is the filtering
    source file for filtering_source, already downloaded and staged by
    Nextflow."""
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
        df = df[["id", context]].rename(columns={context: "expression"})
        return df, None

    if filtering_source == "GTEx":
        df = pd.read_csv(filtering_file, sep="\t", skiprows=2, header=0)
        df.rename(columns={df.columns[0]: "id"}, inplace=True)
        df = df[["id", context]].rename(columns={context: "expression"})
        df["id"] = df["id"].apply(lambda x: x.split(".")[0])
        return df, f"https://gtexportal.org/home/tissue/{context}"
    elif filtering_source == "PAXDB":
        # the first line of the file holds the PaxDB dataset id, e.g. "#id: 4032974247"
        first_line = Path(filtering_file).read_text().splitlines()[0]
        dataset_id = first_line.split(":", 1)[1].strip()
        # the column header itself is also '#'-prefixed in PaxDB files, so it
        # gets dropped as a comment too; supply the column names explicitly
        df = pd.read_csv(
            filtering_file,
            sep="\t",
            comment="#",
            header=None,
            names=["gene_name", "string_external_id", "abundance"],
        )
        df = df.rename(columns={"gene_name": "id", "abundance": "expression"})
        df = df[["id", "expression"]]
        return df, f"https://pax-db.org/dataset/9606/{dataset_id}/"
    elif filtering_source == "TCGA":
        if not context.upper().startswith("TCGA_"):
            raise ValueError(
                f"--context must be of the form 'TCGA_{{cancer_type}}' when using --filtering_source TCGA, got '{context}'"
            )
        df = pd.read_csv(filtering_file, sep="\t", header=0)
        df.rename(columns={df.columns[0]: "id"}, inplace=True)
        df["id"] = df["id"].apply(lambda x: x.split(".")[0])
        sample_columns = df.columns.drop("id")
        # values are log2(tpm+1) transformed; undo that before taking the median
        # so it is computed on the same TPM scale as GTEx
        tpm = 2 ** df[sample_columns] - 1
        df["expression"] = tpm.median(axis=1)
        df = df[["id", "expression"]]
        return df, None
    elif filtering_source == "CRAPome":
        # pass engine explicitly since the file staged by Nextflow may not
        # keep the .xlsx extension (the source URL has a query string)
        raw = pd.read_excel(filtering_file, sheet_name="Sheet1", engine="openpyxl")
        # CC* columns hold the per-experiment spectral counts observed in
        # control (contaminant-background) AP-MS runs; a protein is
        # considered "detected" in an experiment if its spectral count there
        # is greater than zero. The filtering criterion is the fraction of
        # control experiments it was detected in (0-1), not the total
        # spectral count, so the threshold reads as a percentage of experiments
        cc_columns = [col for col in raw.columns if col.startswith("CC")]
        raw["expression"] = (raw[cc_columns].fillna(0) > 0).sum(axis=1) / len(
            cc_columns
        )
        df = raw[["geneSymbol", "expression"]].groupby(
            "geneSymbol", as_index=False
        ).agg({"expression": "sum"})
        df = df.rename(columns={"geneSymbol": "id"})
        # convert to ensembl gene ids here so the returned "id" column matches
        # the same contract as the other sources (native ensembl gene ids,
        # further converted downstream by convert_id_space if needed)
        df = convert_id_space(df, "ensembl", source_space="symbol")
        return df, "https://reprint-apms.org/"
    elif filtering_source == "ProteomicsDB":
        raise NotImplementedError(f"{filtering_source} support is not yet implemented")
    else:
        raise ValueError(
            f"Unknown filtering_source: '{filtering_source}'. Must be one of: GTEx, ProteomicsDB, PAXDB, TCGA, CRAPome"
        )


ID_SPACE_TARGET_NAMESPACES = {
    "ensembl": "ENSG",
    "entrez": "ENTREZGENE_ACC",
    "uniprot": "UNIPROT_GN_ACC",
    "symbol": "HGNC",
}


def convert_id_space(expression_df, id_space, source_space="ensembl"):
    """Converts the "id" column of expression_df, which is assumed to hold
    identifiers in source_space, to id_space. If both spaces are the same, no
    conversion is needed and expression_df is returned unchanged."""
    if id_space == source_space:
        return expression_df
    ids = expression_df["id"]
    if id_space not in ID_SPACE_TARGET_NAMESPACES:
        raise ValueError(f"Unsupported id_space: {id_space}")
    targetSpace = ID_SPACE_TARGET_NAMESPACES[id_space]
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


def save_expression_distribution(
    expression_by_context, stem, context, source, threshold, threshold_label
):
    values = expression_by_context["expression"]
    if values.empty:
        return
    percentiles = list(range(0, 101))
    pct_values = [float(values.quantile(p / 100)) for p in percentiles]
    distribution = {
        "name": f"{stem}.{context}.{source}.{threshold_label}",
        "context": context,
        "source": source,
        "threshold": threshold,
        "data": [[v, p] for p, v in zip(percentiles, pct_values)],
    }
    with open(
        f"{stem}.{context}.{source}.{threshold_label}.expression_distribution.yaml", "w"
    ) as f:
        yaml.safe_dump(distribution, f, sort_keys=False, default_flow_style=None)


def write_crapome_filtering_statistics(
    network_name,
    source,
    threshold_label,
    network_num_vertices,
    nodes_marked_as_crapome,
    genes_filtered_by_threshold,
):
    with open("filtering_statistic.tsv", "w") as f:
        f.write(
            "network\t"
            "nodes_marked_as_crapome_absolute\tnodes_marked_as_crapome_relative\t"
            "genes_filtered_by_threshold_absolute\tgenes_filtered_by_threshold_relative\n"
        )
        f.write(
            f"{network_name}.{source}.{threshold_label}\t"
            f"{nodes_marked_as_crapome}\t"
            f"{nodes_marked_as_crapome / network_num_vertices if network_num_vertices else 0.0}\t"
            f"{genes_filtered_by_threshold}\t"
            f"{genes_filtered_by_threshold / network_num_vertices if network_num_vertices else 0.0}\n"
        )


def write_filtering_statistics(
    network_name,
    context,
    source,
    threshold_label,
    dataset_link,
    n_genes_in_context,
    expression_num_entries,
    network_num_vertices,
    genes_not_in_network,
    genes_not_in_expression_file,
    genes_filtered_by_threshold,
):
    context_label = f"<a href={dataset_link}>{context}</a>" if dataset_link else context
    genes_not_in_network_relative = (
        genes_not_in_network / expression_num_entries if expression_num_entries else 0.0
    )
    genes_not_in_expression_file_relative = (
        genes_not_in_expression_file / network_num_vertices
        if network_num_vertices
        else 0.0
    )
    genes_filtered_by_threshold_relative = (
        genes_filtered_by_threshold / network_num_vertices if network_num_vertices else 0.0
    )
    with open("filtering_statistic.tsv", "w") as f:
        f.write(
            "network\tcontext\tgenes_in_context\t"
            "genes_not_in_network_absolute\tgenes_not_in_network_relative\t"
            "genes_not_in_expression_file_absolute\tgenes_not_in_expression_file_relative\t"
            "genes_filtered_by_threshold_absolute\tgenes_filtered_by_threshold_relative\n"
        )
        f.write(
            f"{network_name}.{context}.{source}.{threshold_label}\t{context_label}\t{n_genes_in_context}\t"
            f"{genes_not_in_network}\t{genes_not_in_network_relative}\t"
            f"{genes_not_in_expression_file}\t{genes_not_in_expression_file_relative}\t"
            f"{genes_filtered_by_threshold}\t{genes_filtered_by_threshold_relative}\n"
        )


def filter_network(
    network_file, threshold_label, expression_by_context, context, source, dataset_link
):
    threshold = float(threshold_label)
    network = gt.load_graph(network_file)
    stem = Path(network_file).stem
    name_index = utils.name2index(network)
    network_num_vertices = network.num_vertices()
    expression_num_entries = expression_by_context.shape[0]
    in_network_expression = expression_by_context[
        expression_by_context["id"].isin(name_index.keys())
    ].copy()
    genes_not_in_network = expression_num_entries - in_network_expression.shape[0]
    genes_not_in_expression_file = (
        network_num_vertices - in_network_expression["id"].nunique()
    )
    # CRAPome scores are contaminant/false-positive likelihoods rather than
    # expression levels, so proteins are kept below the threshold instead of
    # above it.
    is_crapome = source == "CRAPome"
    if is_crapome:
        keep_mask = in_network_expression["expression"] < threshold
    else:
        keep_mask = in_network_expression["expression"] > threshold
    context_specific_in_network = in_network_expression[keep_mask].copy()
    genes_filtered_by_threshold = (
        in_network_expression.shape[0] - context_specific_in_network.shape[0]
    )
    context_specific_in_network["vertex_id"] = context_specific_in_network["id"].map(
        name_index
    )
    in_network_expression["vertex_id"] = in_network_expression["id"].map(name_index)

    # filter network by context specific expression
    network.vp["context_filter"] = network.new_vertex_property("bool")
    network.vp["expression_in_context"] = network.new_vertex_property("double")
    if is_crapome:
        # Contaminant sources (CRAPome) only list a small number of proteins.
        # Absence from the source means "not a known contaminant", so those
        # nodes must be kept. Only the in-network nodes that fail the threshold
        # (the actual contaminants) should be removed; everything else stays.
        network.vp["context_filter"].a = True
        contaminant_vertex_ids = (
            set(in_network_expression["vertex_id"].unique())
            - set(context_specific_in_network["vertex_id"].unique())
        )
        for vertex_id in contaminant_vertex_ids:
            network.vp["context_filter"][vertex_id] = False
        # record the contaminant fraction on the nodes that carry a value
        for vertex_id in in_network_expression["vertex_id"].unique():
            expression_value = in_network_expression.loc[
                in_network_expression["vertex_id"] == vertex_id
            ]["expression"].values[0]
            network.vp["expression_in_context"][vertex_id] = expression_value
    else:
        for vertex_id in context_specific_in_network["vertex_id"].unique():
            network.vp["context_filter"][vertex_id] = True
            expression_value = context_specific_in_network.loc[
                context_specific_in_network["vertex_id"] == vertex_id
            ]["expression"].values[0]
            network.vp["expression_in_context"][vertex_id] = expression_value
    network.set_vertex_filter(network.vp["context_filter"])
    network.purge_vertices()
    network.clear_filters()
    del network.vp["context_filter"]
    
    if is_crapome:
        n_genes_in_context = (expression_by_context["expression"] < threshold).sum()
    else:
        n_genes_in_context = (expression_by_context["expression"] > threshold).sum()
    if str(source).lower() == "crapome":
        network.save(f"{stem}.{source}.{threshold_label}.gt")
        write_crapome_filtering_statistics(
            stem,
            source,
            threshold_label,
            network_num_vertices,
            in_network_expression.shape[0],
            genes_filtered_by_threshold,
        )
        return
    
    save_expression_distribution(
        in_network_expression, stem, context, source, threshold, threshold_label
    )
    network.save(f"{stem}.{context}.{source}.{threshold_label}.gt")
    write_filtering_statistics(
        stem,
        context,
        source,
        threshold_label,
        dataset_link,
        n_genes_in_context,
        expression_num_entries,
        network_num_vertices,
        genes_not_in_network,
        genes_not_in_expression_file,
        genes_filtered_by_threshold,
    )


def main(argv=None):
    args = parse_args(argv)
    global id_space
    id_space = args.id_space
    expression_by_context, dataset_link = resolve_expression_file(
        args.filtering_source,
        args.custom_expression_file,
        args.filtering_file,
        args.context,
    )
    # map gene_ids to the specified id space
    expression_by_context = convert_id_space(expression_by_context, args.id_space)
    # filter network by context specific expression with the given threshold
    filter_network(
        args.network,
        args.threshold,
        expression_by_context,
        args.context,
        args.filtering_source,
        dataset_link,
    )


if __name__ == "__main__":
    sys.exit(main())
