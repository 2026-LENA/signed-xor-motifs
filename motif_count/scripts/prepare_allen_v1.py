#!/usr/bin/env python3
"""Flatten an Allen/Billeh V1 SONATA network for signed_xor_motif_count.c."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as error:  # pragma: no cover - user-facing dependency check
    raise SystemExit(
        "This converter requires numpy and h5py. Install them with "
        "`python3 -m pip install numpy h5py`."
    ) from error


def read_space_table(path: Path) -> list[dict[str, str]]:
    """Read the space-delimited SONATA type tables used by this model."""
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=" ", skipinitialspace=True))


def only_child(group: h5py.Group, description: str) -> h5py.Group:
    names = list(group.keys())
    if len(names) != 1:
        raise RuntimeError(
            f"Expected one {description}, found {len(names)}: {', '.join(names)}"
        )
    child = group[names[0]]
    if not isinstance(child, h5py.Group):
        raise RuntimeError(f"Expected {child.name} to be an HDF5 group")
    return child


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nodes_h5", type=Path)
    parser.add_argument("node_types_csv", type=Path)
    parser.add_argument("edges_h5", type=Path)
    parser.add_argument("edge_types_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--chunk-size", type=int, default=250_000)
    args = parser.parse_args()

    node_type_rows = read_space_table(args.node_types_csv)
    node_type_metadata = {
        int(row["node_type_id"]): (row["pop_name"], row["ei"].lower())
        for row in node_type_rows
    }
    edge_type_rows = read_space_table(args.edge_types_csv)
    max_edge_type = max(int(row["edge_type_id"]) for row in edge_type_rows)
    weight_by_edge_type = np.full(max_edge_type + 1, np.nan, dtype=np.float64)
    for row in edge_type_rows:
        weight_by_edge_type[int(row["edge_type_id"])] = float(row["syn_weight"])

    with h5py.File(args.nodes_h5, "r") as nodes_file:
        node_population = only_child(nodes_file["nodes"], "node population")
        node_ids = node_population["node_id"][:]
        node_type_ids = node_population["node_type_id"][:]

    labels: list[str] = []
    connection_types: list[str] = []
    metadata_code: dict[tuple[str, str], int] = {}
    code_by_node_type: dict[int, int] = {}
    for node_type_id, (label, ei) in node_type_metadata.items():
        connection_type = "EXC" if ei == "e" else "INH" if ei == "i" else "OTHER"
        key = (label, connection_type)
        if key not in metadata_code:
            metadata_code[key] = len(labels)
            labels.append(label)
            connection_types.append(connection_type)
        code_by_node_type[node_type_id] = metadata_code[key]

    max_node_id = int(node_ids.max())
    code_by_node_id = np.full(max_node_id + 1, -1, dtype=np.int16)
    codes = np.fromiter(
        (code_by_node_type[int(node_type_id)] for node_type_id in node_type_ids),
        dtype=np.int16,
        count=len(node_type_ids),
    )
    code_by_node_id[node_ids] = codes

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.edges_h5, "r") as edges_file, args.output_csv.open(
        "w", encoding="utf-8", newline="", buffering=8 * 1024 * 1024
    ) as output_handle:
        edge_population = only_child(edges_file["edges"], "edge population")
        source_dataset = edge_population["source_node_id"]
        target_dataset = edge_population["target_node_id"]
        edge_type_dataset = edge_population["edge_type_id"]
        row_count = len(source_dataset)
        if not (len(target_dataset) == len(edge_type_dataset) == row_count):
            raise RuntimeError("SONATA edge datasets have inconsistent lengths")

        writer = csv.writer(output_handle, lineterminator="\n")
        writer.writerow(
            [
                "source_node_id",
                "target_node_id",
                "neuron_type",
                "connection_type",
                "syn_weight",
            ]
        )

        for start in range(0, row_count, args.chunk_size):
            end = min(start + args.chunk_size, row_count)
            source_ids = source_dataset[start:end]
            target_ids = target_dataset[start:end]
            edge_type_ids = edge_type_dataset[start:end]
            source_codes = code_by_node_id[source_ids]
            weights = weight_by_edge_type[edge_type_ids]
            if np.any(source_codes < 0):
                raise RuntimeError(f"Unknown source node metadata in rows {start}:{end}")
            if np.any(np.isnan(weights)):
                raise RuntimeError(f"Unknown edge type in rows {start}:{end}")

            writer.writerows(
                (
                    int(source),
                    int(target),
                    labels[int(code)],
                    connection_types[int(code)],
                    repr(float(weight)),
                )
                for source, target, code, weight in zip(
                    source_ids, target_ids, source_codes, weights
                )
            )
            print(
                f"\rWrote {end:,}/{row_count:,} rows "
                f"({100.0 * end / row_count:5.1f}%)",
                end="",
                file=sys.stderr,
                flush=True,
            )

    print(f"\nFinished: {args.output_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
