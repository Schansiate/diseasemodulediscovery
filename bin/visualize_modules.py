#!/usr/bin/env python

"""Takes a gt file and saves it as pdf, png, svg, and html."""


import argparse
from collections import defaultdict
import logging
import sys
from pathlib import Path

import util

import graph_tool.all as gt
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import networkx as nx
import pyintergraph
from pyvis.network import Network

logger = logging.getLogger()


def compute_font_size(n, fmt):
    if fmt == "png":
        C = 36
    else:
        C = 18
    font_size = round(C * (25 / n) ** 0.8)
    return max(1, font_size)


def add_drugs(args, g, node_mapping):
    drug_df = pd.read_csv(args.drugs, sep="\t")
    assert "name" in drug_df.columns

    if "drug_name" not in drug_df.columns:
        logger.warning(
            "The drug file does not contain a 'drug_name' column. Skipping drug addition."
        )
        return set()

    drug_df.dropna(subset=["drug_name"], inplace=True)
    drug_df["name"] = drug_df["name"].astype(str)
    drug_df["drug_name"] = drug_df["drug_name"].astype(str)

    gene2drugs = defaultdict(list)
    drug_set = set(drug_df["drug_name"].unique())

    for _, row in drug_df.iterrows():
        if row["name"] not in gene2drugs:
            gene2drugs[row["name"]] = []
        if pd.notna(row["drug_name"]):
            gene2drugs[row["name"]].append(row["drug_name"])
    for name, drug_list in gene2drugs.items():
        if name in node_mapping:
            protein_vertex = node_mapping[name]
        else:
            logger.warning(f"Node {name} not found in the module")
            continue
        for drug in drug_list:
            if drug not in node_mapping:
                drug_vertex = g.add_vertex()
                g.vp["name"][drug_vertex] = drug
                node_mapping[drug] = drug_vertex
            drug_vertex = node_mapping[drug]
            g.add_edge(protein_vertex, drug_vertex)
    return drug_set


def visualize_expression(args, g, pos):
    """Visualize tissue expression values as a PNG, if present."""
    expression_property = "expression_in_tissue"
    if expression_property not in g.vp:
        logger.info(
            "Skipping expression visualization because the graph has no "
            "'expression_in_tissue' vertex property."
        )
        return

    expr = g.vp[expression_property]
    values = np.array([float(expr[v]) for v in g.vertices()], dtype=float)

    # LogNorm requires strictly positive values; floor non-positive
    # expression (e.g. TPM == 0) to a small positive epsilon so it still
    # maps to the low end of the color scale instead of erroring out.
    positive_values = values[values > 0]
    floor = positive_values.min() / 10 if positive_values.size else 1e-3
    plot_values = np.clip(values, floor, None)

    vmin = plot_values.min()
    vmax = plot_values.max()
    if vmin == vmax:
        vmin /= 10
        vmax *= 10

    norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    cmap = mpl.colormaps["Reds"]

    vertex_fill_color = g.new_vertex_property("vector<double>")
    for v, value in zip(g.vertices(), plot_values):
        vertex_fill_color[v] = cmap(norm(value))

    vertex_shape = g.new_vertex_property("string")
    vertex_aspect = g.new_vertex_property("double")
    for v in g.vertices():
        if "is_seed" in g.vp and g.vp["is_seed"][v]:
            vertex_shape[v] = "square"
            vertex_aspect[v] = 1.6
        else:
            vertex_shape[v] = "circle"
            vertex_aspect[v] = 1.0

    raw_output = Path(f"{args.prefix}.expression.raw.png")
    output = f"{args.prefix}.expression.png"
    gt.graph_draw(
        g,
        pos=pos,
        vertex_fill_color=vertex_fill_color,
        vertex_color="black",
        vertex_shape=vertex_shape,
        vertex_aspect=vertex_aspect,
        vertex_size=18,
        bg_color="white",
        output_size=(800, 800),
        output=str(raw_output),
    )

    img = plt.imread(raw_output)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(img)
    ax.axis("off")

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Expression value (log scale)")

    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    raw_output.unlink(missing_ok=True)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse the modules of different tools.",
        epilog="Example: python save_modules.py -m module1.gt -p 'module1'",
    )
    parser.add_argument(
        "-m",
        "--module",
        help="Path to the module output.",
        type=Path,
    )
    parser.add_argument(
        "-d",
        "--drugs",
        help="Path to the drug file.",
        type=Path,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Prefix to name the output files.",
        type=str,
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)

    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.module.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    logger.debug(f"{args=}")

    # load the module file
    g = util.load_graph(str(args.module))
    node_mapping = util.name2index(g)
    if args.drugs:
        drug_set = add_drugs(args, g, node_mapping)

    # color the seed genes red
    g.vp["color"] = g.new_vertex_property("string")
    g.vp["label_node"] = g.new_vertex_property("string")
    g.vp["shape"] = g.new_vertex_property("string")
    g.vp["text_color"] = g.new_vertex_property("string")
    for v in g.vertices():
        g.vp["label_node"][v] = g.vp["name"][v]
        if args.drugs and g.vp["name"][v] in drug_set:
            g.vp["color"][v] = "#90EE90"  # green for drugs
            g.vp["shape"][v] = "triangle"
            g.vp["text_color"][v] = "black"
        else:
            g.vp["shape"][v] = "circle"
            g.vp["text_color"][v] = "white"
            if g.vertex_properties["is_seed"][v]:
                g.vp["color"][v] = "red"  # red for seed genes
            else:
                g.vp["color"][v] = "blue"  # blue for added genes

    # calculate the layout
    pos = gt.sfdp_layout(g)

    visualize_expression(args, g, pos)

    # save as pdf, png, svg
    n = g.num_vertices()
    for format in ["pdf", "png", "svg"]:
        size = (3000, 3000) if format == "png" else (1000, 1000)
        font_size = compute_font_size(n, format)
        gt.graph_draw(
            g,
            pos,
            vertex_fill_color=g.vp["color"],
            vertex_text_color=g.vp["text_color"],
            output_size=size,
            vertex_text=g.vp["label_node"],
            vorder=g.vp["is_seed"],
            vertex_shape=g.vp["shape"],
            vertex_font_size=font_size,
            edge_pen_width=0.5,
            output=f"{args.prefix}.{format}",
            vertex_size=1,
            vertex_pen_width=0.5,
            edge_color="#A9A9A9",
        )

    # add x and y properties for consistent html layout
    g.vp["x"] = g.new_vertex_property("double")
    g.vp["y"] = g.new_vertex_property("double")

    for v in g.vertices():
        g.vp["x"][v] = pos[v][0] * 30
        g.vp["y"][v] = pos[v][1] * 30

    # convert to networkx and save as html with pyvis
    nx_graph = pyintergraph.gt2nx(g, labelname="label_node")

    nt = Network()
    nt.from_nx(nx_graph)

    # get vertex properties as dataframe
    vp_df = util.vp2df(g)

    # add titles
    for node in nt.nodes:
        row = vp_df.loc[node["name"]]
        node_properties = "\n".join(f"{col}: {val}" for col, val in row.items())
        node[
            "title"
        ] = f"""
            {node['name']}\n
            {node_properties}
        """

    # turn off physics
    nt.toggle_physics(False)

    # show setting buttons
    nt.show_buttons(filter_=["physics", "interaction", "manipulation"])
    # save as html
    nt.show(f"{args.prefix}.html")


if __name__ == "__main__":
    sys.exit(main())
