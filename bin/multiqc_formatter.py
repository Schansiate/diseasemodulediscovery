#!/usr/bin/env python
import argparse
from pathlib import Path
import yaml
import sys


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="formats file for multiqc custom contents",
        epilog="Example: python multiqc_formatter.py -i network.gt -f network_degree",
    )
    parser.add_argument(
        "-i", "--input", type=Path, nargs="*", required=True, help="Input files"
    )
    parser.add_argument("-H", "--header", type=Path, required=True, help="Header file")
    return parser.parse_args(argv)


def parse_input(input_files, header_file):
    with open(header_file, "r") as header:
        header_data = yaml.safe_load(header)
        header_id = header_data.get("id", "")

    if header_id == "network_node_degree_distribution":
        save_node_degree_distribution(input_files, header_file)
    elif header_id == "filtered_network_expression_distribution":
        save_expression_distribution(input_files, header_file)


def save_node_degree_distribution(input_files, header_file):
    with open(header_file, "r", encoding="utf-8") as header:
        mqc_payload = yaml.safe_load(header) or {}

    absolute_data = {}
    relative_data = {}

    for file in input_files:
        with open(file, "r", encoding="utf-8") as distribution_file:
            distribution = yaml.safe_load(distribution_file) or {}

        network_name = distribution.get("name") or file.stem
        absolute = distribution.get("absolute")
        relative = distribution.get("relative")

        if absolute is None or relative is None:
            raise ValueError(
                f"Invalid distribution YAML in {file}: expected keys 'absolute' and 'relative'"
            )

        absolute_data[network_name] = absolute
        relative_data[network_name] = relative

    mqc_payload["data"] = [absolute_data, relative_data]

    with open("./node_degree_distribution_mqc.yaml", "w", encoding="utf-8") as file:
        yaml.safe_dump(mqc_payload, file, default_flow_style=None)


def save_expression_distribution(input_files, header_file):
    with open(header_file, "r", encoding="utf-8") as header:
        mqc_payload = yaml.safe_load(header) or {}

    data = {}
    thresholds_by_value = {}

    for file in input_files:
        with open(file, "r", encoding="utf-8") as distribution_file:
            distribution = yaml.safe_load(distribution_file) or {}

        dist_data = distribution.get("data")

        if dist_data is None:
            raise ValueError(
                f"Invalid distribution YAML in {file}: expected key 'data'"
            )

        context = distribution.get("context")
        source = distribution.get("source")
        threshold = distribution.get("threshold")

        if context is not None and source is not None and threshold is not None:
            line_name = f"{context}.{source}.{threshold}"
        else:
            line_name = distribution.get("name") or file.stem

        data[line_name] = dist_data

        if threshold is not None:
            thresholds_by_value.setdefault(threshold, []).append(line_name)

    if thresholds_by_value:
        x_lines = []
        for value, names in sorted(thresholds_by_value.items()):
            label = f"Threshold ({value})" if len(names) == len(input_files) else (
                f"Threshold ({value}): {', '.join(names)}"
            )
            x_lines.append(
                {
                    "value": value,
                    "width": 2,
                    "dash": "dash",
                    "label": label,
                }
            )
        mqc_payload.setdefault("pconfig", {})["x_lines"] = x_lines

    mqc_payload["data"] = data

    with open("./expression_distribution_mqc.yaml", "w", encoding="utf-8") as file:
        yaml.safe_dump(mqc_payload, file, sort_keys=False, default_flow_style=None)


def main():
    args = parse_args()
    parse_input(args.input, args.header)


if __name__ == "__main__":
    sys.exit(main())
