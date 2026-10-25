#!/usr/bin/env python3
"""Audit FlyWire annotation coverage and class-block swap viability."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile


DEFAULT_FIELDS = ("super_class", "cell_class", "cell_sub_class", "cell_type")
MISSING_LABEL = "__UNANNOTATED__"


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
            stream.write(text)
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def canonical_sign(nt_type: str) -> str | None:
    value = nt_type.strip().upper()
    if value == "ACH":
        return "EXC"
    if value in {"GABA", "GLUT"}:
        return "INH"
    return None


def read_annotations(
    path: Path, fields: tuple[str, ...]
) -> tuple[dict[str, tuple[str, ...]], dict[str, int], int]:
    annotations: dict[str, tuple[str, ...]] = {}
    populated = {field: 0 for field in fields}
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        available = set(reader.fieldnames or ())
        required = {"root_id", *fields}
        missing = sorted(required - available)
        if missing:
            raise ValueError(
                f"annotation table is missing columns: {', '.join(missing)}"
            )
        rows = 0
        for row_number, row in enumerate(reader, start=2):
            rows += 1
            root_id = (row.get("root_id") or "").strip()
            if not root_id:
                raise ValueError(f"empty root_id at annotation row {row_number}")
            if root_id in annotations:
                raise ValueError(f"duplicate root_id {root_id} in annotation table")
            labels = tuple(sys.intern((row.get(field) or "").strip()) for field in fields)
            annotations[root_id] = labels
            for field, label in zip(fields, labels):
                if label:
                    populated[field] += 1
    return annotations, populated, rows


def graph_nodes_and_counts(path: Path) -> tuple[dict[str, int], dict[str, int]]:
    node_ids: dict[str, int] = {}
    counts = {"rows": 0, "signed_edges": 0, "other_edges": 0, "self_loops": 0}
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"pre_root_id", "post_root_id", "nt_type"}
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"graph is missing columns: {', '.join(missing)}")
        for row_number, row in enumerate(reader, start=2):
            pre = (row.get("pre_root_id") or "").strip()
            post = (row.get("post_root_id") or "").strip()
            nt_type = row.get("nt_type") or ""
            if not pre or not post or not nt_type.strip():
                raise ValueError(f"empty pre/post/nt_type at graph row {row_number}")
            if pre not in node_ids:
                node_ids[pre] = len(node_ids)
            if post not in node_ids:
                node_ids[post] = len(node_ids)
            counts["rows"] += 1
            if pre == post:
                counts["self_loops"] += 1
            if canonical_sign(nt_type) is None:
                counts["other_edges"] += 1
            else:
                counts["signed_edges"] += 1
    counts["nodes"] = len(node_ids)
    return node_ids, counts


def audit_field(
    graph: Path,
    annotations: dict[str, tuple[str, ...]],
    node_ids: dict[str, int],
    field_index: int,
    field: str,
    annotation_populated_rows: int,
    annotation_rows: int,
) -> dict:
    labels: dict[str, str] = {}
    graph_labels: set[str] = set()
    annotated_nodes = 0
    for root_id in node_ids:
        record = annotations.get(root_id)
        label = record[field_index] if record else ""
        if label:
            annotated_nodes += 1
            graph_labels.add(label)
            labels[root_id] = label

    blocks: dict[tuple[str, str, str], list] = {}
    both_annotated_edges = 0
    with graph.open("r", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            sign = canonical_sign(row["nt_type"])
            if sign is None:
                continue
            pre = row["pre_root_id"].strip()
            post = row["post_root_id"].strip()
            pre_label = labels.get(pre, MISSING_LABEL)
            post_label = labels.get(post, MISSING_LABEL)
            if pre_label != MISSING_LABEL and post_label != MISSING_LABEL:
                both_annotated_edges += 1
            key = (sign, pre_label, post_label)
            block = blocks.get(key)
            if block is None:
                block = [0, set(), set()]
                blocks[key] = block
            block[0] += 1
            block[1].add(node_ids[pre])
            block[2].add(node_ids[post])

    edge_total = sum(block[0] for block in blocks.values())
    eligible_blocks = 0
    eligible_edges = 0
    singleton_blocks = 0
    singleton_edges = 0
    blocks_ge_10 = 0
    edges_ge_10 = 0
    blocks_ge_100 = 0
    edges_ge_100 = 0
    largest_block = 0
    for edge_count, sources, targets in blocks.values():
        largest_block = max(largest_block, edge_count)
        if edge_count == 1:
            singleton_blocks += 1
            singleton_edges += 1
        if edge_count >= 10:
            blocks_ge_10 += 1
            edges_ge_10 += edge_count
        if edge_count >= 100:
            blocks_ge_100 += 1
            edges_ge_100 += edge_count
        if edge_count >= 2 and len(sources) >= 2 and len(targets) >= 2:
            eligible_blocks += 1
            eligible_edges += edge_count

    edge_denominator = max(edge_total, 1)
    nodes = max(len(node_ids), 1)
    annotation_denominator = max(annotation_rows, 1)
    return {
        "field": field,
        "annotation_rows_with_value": annotation_populated_rows,
        "annotation_row_coverage": annotation_populated_rows
        / annotation_denominator,
        "graph_nodes_with_value": annotated_nodes,
        "graph_node_coverage": annotated_nodes / nodes,
        "distinct_graph_labels": len(graph_labels),
        "signed_edges_with_both_endpoints_annotated": both_annotated_edges,
        "signed_edge_annotation_coverage": both_annotated_edges
        / edge_denominator,
        "blocks": len(blocks),
        "eligible_blocks": eligible_blocks,
        "eligible_edges": eligible_edges,
        "eligible_edge_fraction": eligible_edges / edge_denominator,
        "frozen_edges": edge_total - eligible_edges,
        "singleton_blocks": singleton_blocks,
        "singleton_edges": singleton_edges,
        "blocks_with_at_least_10_edges": blocks_ge_10,
        "edges_in_blocks_with_at_least_10_edges": edges_ge_10,
        "blocks_with_at_least_100_edges": blocks_ge_100,
        "edges_in_blocks_with_at_least_100_edges": edges_ge_100,
        "largest_block_edges": largest_block,
    }


def percent(value: float) -> str:
    return f"{100.0 * value:.2f}%"


def render_markdown(report: dict) -> str:
    rows = []
    for result in report["fields"]:
        rows.append(
            "| {field} | {labels:,} | {nodes} | {edges} | {blocks:,} | "
            "{eligible} | {frozen:,} | {largest:,} |".format(
                field=result["field"],
                labels=result["distinct_graph_labels"],
                nodes=percent(result["graph_node_coverage"]),
                edges=percent(result["signed_edge_annotation_coverage"]),
                blocks=result["blocks"],
                eligible=percent(result["eligible_edge_fraction"]),
                frozen=result["frozen_edges"],
                largest=result["largest_block_edges"],
            )
        )
    table = "\n".join(rows)
    chosen = report["recommended_field"] or "none"
    return f"""# FlyWire v783 class-block audit

| Input | SHA-256 |
|---|---|
| `{report['graph']['path']}` | `{report['graph']['sha256']}` |
| `{report['annotations']['path']}` | `{report['annotations']['sha256']}` |

The graph contains {report['graph']['nodes']:,} neurons,
{report['graph']['signed_edges']:,} ACH/GABA/GLUT signed edges, and
{report['graph']['other_edges']:,} fixed other-transmitter edges. The annotation
table contains {report['annotations']['rows']:,} neurons.

| Field | Labels in graph | Node coverage | Signed-edge coverage | Sign/class blocks | Structurally eligible edges | Frozen edges | Largest block |
|---|---:|---:|---:|---:|---:|---:|---:|
{table}

"Structurally eligible" means an edge lies in a sign/source-class/target-class
block with at least two edges, two distinct sources, and two distinct targets.
It is an upper bound on edges that will actually move because duplicate-edge
and self-loop rejection can still prevent individual swaps.

Recommended class field: **`{chosen}`**.

The recommendation selects the first (coarsest) requested hierarchy level with
at least {percent(report['selection_thresholds']['minimum_node_coverage'])}
graph-node coverage and {percent(report['selection_thresholds']['minimum_eligible_edge_fraction'])}
structurally eligible signed edges. Unannotated neurons are assigned to a
single `{MISSING_LABEL}` block so signed in/out degree remains preserved for
the full graph.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graph", type=Path, help="sign-normalized FlyWire CSV")
    parser.add_argument("annotations", type=Path, help="FlyWire annotation TSV")
    parser.add_argument("--fields", nargs="+", default=list(DEFAULT_FIELDS))
    parser.add_argument("--report-json", type=Path)
    parser.add_argument("--report-md", type=Path)
    parser.add_argument("--minimum-node-coverage", type=float, default=0.9)
    parser.add_argument("--minimum-eligible-edge-fraction", type=float, default=0.95)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    fields = tuple(args.fields)
    if len(set(fields)) != len(fields):
        print("error: --fields contains duplicates", file=sys.stderr)
        return 2
    for name, value in (
        ("--minimum-node-coverage", args.minimum_node_coverage),
        ("--minimum-eligible-edge-fraction", args.minimum_eligible_edge_fraction),
    ):
        if not 0.0 <= value <= 1.0:
            print(f"error: {name} must be between zero and one", file=sys.stderr)
            return 2
    try:
        graph = args.graph.resolve()
        annotations_path = args.annotations.resolve()
        if not graph.is_file():
            raise FileNotFoundError(graph)
        if not annotations_path.is_file():
            raise FileNotFoundError(annotations_path)
        annotations, populated, annotation_rows = read_annotations(
            annotations_path, fields
        )
        node_ids, graph_counts = graph_nodes_and_counts(graph)
        results = []
        for index, field in enumerate(fields):
            print(f"Auditing {field}...", file=sys.stderr, flush=True)
            results.append(
                audit_field(
                    graph,
                    annotations,
                    node_ids,
                    index,
                    field,
                    populated[field],
                    annotation_rows,
                )
            )
        recommended = next(
            (
                result["field"]
                for result in results
                if result["graph_node_coverage"] >= args.minimum_node_coverage
                and result["eligible_edge_fraction"]
                >= args.minimum_eligible_edge_fraction
            ),
            None,
        )
        report = {
            "schema_version": 1,
            "graph": {
                "path": str(graph),
                "sha256": sha256_file(graph),
                **graph_counts,
            },
            "annotations": {
                "path": str(annotations_path),
                "sha256": sha256_file(annotations_path),
                "rows": annotation_rows,
            },
            "missing_label": MISSING_LABEL,
            "selection_thresholds": {
                "minimum_node_coverage": args.minimum_node_coverage,
                "minimum_eligible_edge_fraction": args.minimum_eligible_edge_fraction,
            },
            "recommended_field": recommended,
            "fields": results,
        }
        markdown = render_markdown(report)
        if args.report_json:
            atomic_text(
                args.report_json.resolve(),
                json.dumps(report, indent=2, sort_keys=True) + "\n",
            )
        if args.report_md:
            atomic_text(args.report_md.resolve(), markdown)
        if not args.report_json and not args.report_md:
            print(markdown, end="")
        return 0
    except (OSError, ValueError, csv.Error) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
