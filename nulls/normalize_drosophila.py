#!/usr/bin/env python3
"""Collapse repeated FlyWire rows into deterministic binary edges.

FlyWire may report the same neuron-to-neuron connection in several neuropils.
The motif counter is a binary topological search, so feeding those rows
directly to it introduces edge multiplicity.  By default, this utility creates
one edge per (pre, post, nt_type).  With --collapse sign, it instead creates
one edge per (pre, post, sign) and writes canonical ACH, GABA, or OTHER types.
The latter is the appropriate input for a strictly binary signed-topology
analysis.

The output is deterministic and disk-backed.  syn_count is summed and the
neuropil field is set to AGGREGATED.

Only the Python standard library is required.  SQLite keeps memory use bounded
for multi-million-row inputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import sys
import tempfile
import time


REQUIRED_COLUMNS = (
    "pre_root_id",
    "post_root_id",
    "neuropil",
    "syn_count",
    "nt_type",
)


def edge_type(nt_type: str, collapse: str) -> str:
    if collapse == "nt":
        return nt_type
    if collapse != "sign":
        raise ValueError(f"unsupported collapse key: {collapse}")
    if nt_type == "ACH":
        return "ACH"
    if nt_type in {"GABA", "GLUT"}:
        return "GABA"
    return "OTHER"


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def normalize(
    input_path: Path,
    output_path: Path,
    report_path: Path,
    *,
    force: bool = False,
    batch_size: int = 20_000,
    collapse: str = "nt",
) -> dict:
    input_path = input_path.resolve()
    output_path = output_path.resolve()
    report_path = report_path.resolve()

    if len({input_path, output_path, report_path}) != 3:
        raise ValueError("input, output, and report paths must be different")
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if output_path.exists() and not force:
        raise FileExistsError(f"output already exists: {output_path} (use --force)")
    if report_path.exists() and not force:
        raise FileExistsError(f"report already exists: {report_path} (use --force)")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if collapse not in {"nt", "sign"}:
        raise ValueError("collapse must be 'nt' or 'sign'")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    db_fd, db_name = tempfile.mkstemp(
        prefix=".drosophila-normalize-", suffix=".sqlite", dir=output_path.parent
    )
    os.close(db_fd)
    out_fd, out_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", suffix=".tmp", dir=output_path.parent
    )
    os.close(out_fd)

    started = time.monotonic()
    original_rows = 0
    connection_synapses = 0
    connection = sqlite3.connect(db_name)
    try:
        connection.executescript(
            """
            PRAGMA journal_mode = OFF;
            PRAGMA synchronous = OFF;
            PRAGMA temp_store = FILE;
            PRAGMA cache_size = -262144;
            CREATE TABLE edges (
                id INTEGER PRIMARY KEY,
                pre_root_id TEXT NOT NULL,
                post_root_id TEXT NOT NULL,
                nt_type TEXT NOT NULL,
                syn_count INTEGER NOT NULL,
                source_rows INTEGER NOT NULL DEFAULT 1,
                UNIQUE (pre_root_id, post_root_id, nt_type)
            );
            """
        )
        insert_sql = """
            INSERT INTO edges
                (pre_root_id, post_root_id, nt_type, syn_count)
            VALUES (?, ?, ?, ?)
            ON CONFLICT (pre_root_id, post_root_id, nt_type) DO UPDATE SET
                syn_count = edges.syn_count + excluded.syn_count,
                source_rows = edges.source_rows + 1
        """

        batch: list[tuple[str, str, str, int]] = []
        connection.execute("BEGIN")
        with input_path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            missing = [name for name in REQUIRED_COLUMNS if name not in (reader.fieldnames or ())]
            if missing:
                raise ValueError(f"missing required CSV columns: {', '.join(missing)}")
            for line_number, row in enumerate(reader, start=2):
                pre = (row["pre_root_id"] or "").strip()
                post = (row["post_root_id"] or "").strip()
                source_nt_type = (row["nt_type"] or "").strip().upper()
                if not pre or not post or not source_nt_type:
                    raise ValueError(f"empty pre/post/nt_type at CSV line {line_number}")
                nt_type = edge_type(source_nt_type, collapse)
                try:
                    syn_count = int((row["syn_count"] or "").strip())
                except ValueError as error:
                    raise ValueError(
                        f"invalid integer syn_count at CSV line {line_number}: "
                        f"{row['syn_count']!r}"
                    ) from error
                if syn_count < 0:
                    raise ValueError(f"negative syn_count at CSV line {line_number}")
                batch.append((pre, post, nt_type, syn_count))
                original_rows += 1
                connection_synapses += syn_count
                if len(batch) >= batch_size:
                    connection.executemany(insert_sql, batch)
                    batch.clear()
                    if original_rows % 1_000_000 < batch_size:
                        print(f"normalized {original_rows:,} input rows...", file=sys.stderr)
            if batch:
                connection.executemany(insert_sql, batch)
        connection.commit()

        unique_edges = connection.execute(
            "SELECT COUNT(*) FROM edges"
        ).fetchone()[0]
        output_synapses = connection.execute(
            "SELECT COALESCE(SUM(syn_count), 0) FROM edges"
        ).fetchone()[0]
        unique_directed_pairs = connection.execute(
            "SELECT COUNT(*) FROM "
            "(SELECT 1 FROM edges GROUP BY pre_root_id, post_root_id)"
        ).fetchone()[0]
        multi_nt_pairs = connection.execute(
            "SELECT COUNT(*) FROM "
            "(SELECT 1 FROM edges GROUP BY pre_root_id, post_root_id "
            " HAVING COUNT(*) > 1)"
        ).fetchone()[0]
        maximum_source_rows = connection.execute(
            "SELECT COALESCE(MAX(source_rows), 0) FROM edges"
        ).fetchone()[0]

        with open(out_name, "w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(REQUIRED_COLUMNS)
            cursor = connection.execute(
                "SELECT pre_root_id, post_root_id, syn_count, nt_type "
                "FROM edges ORDER BY id"
            )
            for pre, post, syn_count, nt_type in cursor:
                writer.writerow((pre, post, "AGGREGATED", syn_count, nt_type))

        if connection_synapses != output_synapses:
            raise RuntimeError(
                "internal error: total syn_count changed during normalization"
            )
        os.replace(out_name, output_path)

        report = {
            "schema_version": 2,
            "operation": f"collapse_by_pre_post_{'sign' if collapse == 'sign' else 'nt_type'}",
            "collapse_key": collapse,
            "input": str(input_path),
            "input_sha256": sha256_file(input_path),
            "output": str(output_path),
            "output_sha256": sha256_file(output_path),
            "original_rows": original_rows,
            "unique_edges": unique_edges,
            "unique_directed_pairs": unique_directed_pairs,
            "collapsed_duplicate_rows": original_rows - unique_edges,
            "directed_pairs_with_multiple_output_types": multi_nt_pairs,
            "maximum_rows_collapsed_into_one_edge": maximum_source_rows,
            "total_syn_count": connection_synapses,
            "neuropil_policy": "AGGREGATED",
            "output_type_policy": "canonical_sign" if collapse == "sign" else "nt_type",
            "elapsed_seconds": time.monotonic() - started,
        }
        atomic_json(report_path, report)
        return report
    finally:
        connection.close()
        for temporary in (db_name, out_name):
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="original FlyWire CSV")
    parser.add_argument("output", type=Path, help="normalized CSV to create")
    parser.add_argument(
        "--report",
        type=Path,
        help="JSON report path (default: OUTPUT.normalization.json)",
    )
    parser.add_argument("--force", action="store_true", help="replace existing output/report")
    parser.add_argument("--batch-size", type=int, default=20_000)
    parser.add_argument(
        "--collapse",
        choices=("nt", "sign"),
        default="nt",
        help="edge identity to preserve (default: nt)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report_path = args.report or Path(f"{args.output}.normalization.json")
    try:
        report = normalize(
            args.input,
            args.output,
            report_path,
            force=args.force,
            batch_size=args.batch_size,
            collapse=args.collapse,
        )
    except (OSError, ValueError, sqlite3.Error, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
