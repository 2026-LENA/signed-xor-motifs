#!/usr/bin/env python3
"""Run and summarize a resumable Drosophila signed-XOR null experiment.

The runner counts the normalized observed graph, generates one randomized CSV
at a time, counts it with the CUDA executable, records a checkpoint, and then
removes the large temporary CSV by default.  Re-running with a larger
``--nulls`` value extends the same seed range without repeating completed work.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import statistics
import subprocess
import sys
import tempfile
import time


COUNT_PATTERN = re.compile(r"^MOTIF_COUNT_RAW\s+(\d+)", re.MULTILINE)
MANIFEST_VERSION = 1


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


def atomic_json(path: Path, value: dict) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def parse_motif_count(output: str) -> int:
    matches = COUNT_PATTERN.findall(output)
    if len(matches) != 1:
        raise ValueError(
            f"expected one MOTIF_COUNT_RAW line, found {len(matches)}"
        )
    return int(matches[0])


def summarize(observed: int, null_counts: list[int]) -> dict:
    if observed < 0 or any(value < 0 for value in null_counts):
        raise ValueError("motif counts cannot be negative")
    n = len(null_counts)
    if not n:
        return {
            "observed": observed,
            "n_nulls": 0,
            "null_mean": None,
            "null_sample_sd": None,
            "fold_enrichment": None,
            "z_score": None,
            "p_enrichment": None,
            "p_depletion": None,
            "p_two_sided": None,
            "minimum_attainable_empirical_p": None,
            "null_min": None,
            "null_max": None,
        }

    mean = statistics.fmean(null_counts)
    sample_sd = statistics.stdev(null_counts) if n >= 2 else None
    fold = observed / mean if mean else math.inf if observed else math.nan
    z_score = (
        (observed - mean) / sample_sd
        if sample_sd is not None and sample_sd > 0
        else None
    )
    p_enrichment = (1 + sum(value >= observed for value in null_counts)) / (n + 1)
    p_depletion = (1 + sum(value <= observed for value in null_counts)) / (n + 1)
    return {
        "observed": observed,
        "n_nulls": n,
        "null_mean": mean,
        "null_sample_sd": sample_sd,
        "fold_enrichment": fold,
        "z_score": z_score,
        "p_enrichment": p_enrichment,
        "p_depletion": p_depletion,
        "p_two_sided": min(1.0, 2.0 * min(p_enrichment, p_depletion)),
        "minimum_attainable_empirical_p": 1.0 / (n + 1),
        "null_min": min(null_counts),
        "null_max": max(null_counts),
    }


def format_value(value: object) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if math.isnan(value):
            return "undefined"
        if math.isinf(value):
            return "infinity"
        return f"{value:.8g}"
    return str(value)


def render_summary(summary: dict, manifest: dict) -> str:
    rows = (
        ("Observed count", summary["observed"]),
        ("Completed nulls", summary["n_nulls"]),
        ("Null mean", summary["null_mean"]),
        ("Null sample SD", summary["null_sample_sd"]),
        ("Fold enrichment", summary["fold_enrichment"]),
        ("Z-score", summary["z_score"]),
        ("One-sided empirical p (enrichment)", summary["p_enrichment"]),
        ("One-sided empirical p (depletion)", summary["p_depletion"]),
        ("Two-sided empirical p", summary["p_two_sided"]),
        ("Smallest attainable empirical p", summary["minimum_attainable_empirical_p"]),
        ("Null minimum", summary["null_min"]),
        ("Null maximum", summary["null_max"]),
    )
    table = "\n".join(f"| {name} | {format_value(value)} |" for name, value in rows)
    induced = "yes" if manifest["induced"] else "no"
    return f"""# Drosophila signed-XOR null experiment

This file is generated after every completed null count.

| Setting | Value |
|---|---|
| Input | `{manifest['input']}` |
| Input SHA-256 | `{manifest['input_sha256']}` |
| Randomization mode | `{manifest['mode']}` |
| Mixing attempts per edge | {manifest['mixing']} |
| Other neurotransmitters randomized | {manifest['include_other']} |
| Induced motif search | {induced} |

| Statistic | Value |
|---|---:|
{table}

The enrichment p-value uses `(1 + number of null counts >= observed) / (N + 1)`.
The Z-score uses the sample standard deviation of the null counts.  Treat the
statistics as preliminary until the requested number of nulls has completed.
"""


def run_logged(command: list[str], log_path: Path, env: dict[str, str]) -> str:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command_line = shlex.join(command)
    chunks: list[str] = []
    with log_path.open("w", encoding="utf-8", newline="") as log:
        heading = f"$ {command_line}\n"
        print(heading, end="", flush=True)
        log.write(heading)
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=env,
        )
        assert process.stdout is not None
        try:
            for line in process.stdout:
                print(line, end="", flush=True)
                log.write(line)
                chunks.append(line)
            return_code = process.wait()
        except KeyboardInterrupt:
            process.send_signal(signal.SIGINT)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.terminate()
            raise
    if return_code:
        raise RuntimeError(
            f"command failed with exit code {return_code}; see {log_path}"
        )
    return "".join(chunks)


def read_counts(path: Path) -> dict[int, dict]:
    if not path.exists():
        return {}
    rows: dict[int, dict] = {}
    with path.open("r", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            seed = int(row["seed"])
            if seed in rows:
                raise ValueError(f"duplicate seed {seed} in {path}")
            rows[seed] = {
                "seed": seed,
                "count": int(row["count"]),
                "randomize_seconds": float(row["randomize_seconds"]),
                "count_seconds": float(row["count_seconds"]),
            }
    return rows


def write_counts(path: Path, rows: dict[int, dict]) -> None:
    lines = ["seed,count,randomize_seconds,count_seconds\n"]
    for seed in sorted(rows):
        row = rows[seed]
        lines.append(
            f"{seed},{row['count']},{row['randomize_seconds']:.6f},"
            f"{row['count_seconds']:.6f}\n"
        )
    atomic_text(path, "".join(lines))


def executable_identity(path: Path) -> dict:
    resolved = path.resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    if not os.access(resolved, os.X_OK):
        raise PermissionError(f"not executable: {resolved}")
    return {"path": str(resolved), "sha256": sha256_file(resolved)}


def build_manifest(args: argparse.Namespace) -> dict:
    input_path = args.input.resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    return {
        "schema_version": MANIFEST_VERSION,
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "randomizer": executable_identity(args.randomizer),
        "counter": executable_identity(args.counter),
        "mode": args.mode,
        "mixing": args.mixing,
        "include_other": args.include_other,
        "induced": args.induced,
        "seed_start": args.seed_start,
    }


def load_or_create_manifest(path: Path, expected: dict) -> dict:
    if path.exists():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != expected:
            differences = [
                key for key in sorted(set(existing) | set(expected))
                if existing.get(key) != expected.get(key)
            ]
            raise ValueError(
                "experiment configuration differs from its manifest "
                f"({', '.join(differences)}); use a new --output-dir"
            )
        return existing
    atomic_json(path, expected)
    return expected


def update_summary(output_dir: Path, observed: int, rows: dict[int, dict], manifest: dict) -> None:
    values = [rows[seed]["count"] for seed in sorted(rows)]
    result = summarize(observed, values)
    atomic_json(output_dir / "summary.json", result)
    atomic_text(output_dir / "summary.md", render_summary(result, manifest))


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="normalized FlyWire CSV")
    parser.add_argument("--nulls", type=int, required=True, help="total nulls desired")
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--mixing", type=int, default=100)
    parser.add_argument("--mode", choices=("by-nt", "lump-inhibitory"), default="by-nt")
    parser.add_argument("--include-other", action="store_true")
    parser.add_argument("--induced", action="store_true")
    parser.add_argument(
        "--randomizer",
        type=Path,
        default=root / "nulls" / "randomize_drosophila",
    )
    parser.add_argument(
        "--counter",
        type=Path,
        default=root / "motif_count" / "build" / "cuda" / "drosophila_signed_xor_motif",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "motif_count" / "results" / "drosophila_nulls",
    )
    parser.add_argument("--device", help="CUDA device index for CUDA_VISIBLE_DEVICES")
    parser.add_argument("--keep-null-csv", action="store_true")
    parser.add_argument(
        "--observed-count",
        type=int,
        help="trusted normalized observed count; otherwise count input automatically",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.nulls < 1:
        print("error: --nulls must be positive", file=sys.stderr)
        return 2
    if args.seed_start < 0:
        print("error: --seed-start cannot be negative", file=sys.stderr)
        return 2
    if args.mixing < 1:
        print("error: --mixing must be positive", file=sys.stderr)
        return 2
    if args.observed_count is not None and args.observed_count < 0:
        print("error: --observed-count cannot be negative", file=sys.stderr)
        return 2

    output_dir = args.output_dir.resolve()
    logs_dir = output_dir / "logs"
    work_dir = output_dir / "work"
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    if args.device is not None:
        env["CUDA_VISIBLE_DEVICES"] = args.device

    try:
        expected_manifest = build_manifest(args)
        manifest = load_or_create_manifest(output_dir / "manifest.json", expected_manifest)
        observed_path = output_dir / "observed.json"
        counter = manifest["counter"]["path"]
        randomizer = manifest["randomizer"]["path"]
        input_path = manifest["input"]
        counter_options = ["--gpu"] + (["--induced"] if args.induced else [])

        if observed_path.exists():
            observed_record = json.loads(observed_path.read_text(encoding="utf-8"))
            observed = int(observed_record["count"])
            if args.observed_count is not None and args.observed_count != observed:
                raise ValueError("--observed-count differs from the saved observed count")
            print(f"Observed count already complete: {observed:,}")
        elif args.observed_count is not None:
            observed = args.observed_count
            atomic_json(
                observed_path,
                {"count": observed, "source": "provided", "input_sha256": manifest["input_sha256"]},
            )
            print(f"Recorded provided observed count: {observed:,}")
        else:
            started = time.monotonic()
            output = run_logged(
                [counter, *counter_options, input_path],
                logs_dir / "observed.log",
                env,
            )
            observed = parse_motif_count(output)
            atomic_json(
                observed_path,
                {
                    "count": observed,
                    "source": "counted",
                    "input_sha256": manifest["input_sha256"],
                    "elapsed_seconds": time.monotonic() - started,
                },
            )

        counts_path = output_dir / "null_counts.csv"
        rows = read_counts(counts_path)
        update_summary(output_dir, observed, rows, manifest)
        seeds = range(args.seed_start, args.seed_start + args.nulls)
        for index, seed in enumerate(seeds, start=1):
            if seed in rows:
                print(f"[{index}/{args.nulls}] seed {seed}: already complete")
                continue

            print(f"[{index}/{args.nulls}] seed {seed}: randomizing")
            final_csv = work_dir / f"null_{seed}.csv"
            partial_csv = work_dir / f"null_{seed}.csv.partial"
            partial_csv.unlink(missing_ok=True)
            if final_csv.exists() and not args.keep_null_csv:
                final_csv.unlink()
            randomizer_command = [
                randomizer,
                input_path,
                "--out",
                str(partial_csv),
                "--seed",
                str(seed),
                "--mixing",
                str(args.mixing),
                f"--{args.mode}",
            ]
            if args.include_other:
                randomizer_command.append("--include-other")
            randomize_started = time.monotonic()
            run_logged(
                randomizer_command,
                logs_dir / f"null_{seed}_randomize.log",
                env,
            )
            randomize_seconds = time.monotonic() - randomize_started
            os.replace(partial_csv, final_csv)

            print(f"[{index}/{args.nulls}] seed {seed}: counting")
            count_started = time.monotonic()
            output = run_logged(
                [counter, *counter_options, str(final_csv)],
                logs_dir / f"null_{seed}_count.log",
                env,
            )
            count_seconds = time.monotonic() - count_started
            count = parse_motif_count(output)
            rows[seed] = {
                "seed": seed,
                "count": count,
                "randomize_seconds": randomize_seconds,
                "count_seconds": count_seconds,
            }
            write_counts(counts_path, rows)
            update_summary(output_dir, observed, rows, manifest)
            if not args.keep_null_csv:
                final_csv.unlink()
            print(f"[{index}/{args.nulls}] seed {seed}: count={count:,} complete")

        update_summary(output_dir, observed, rows, manifest)
        print(f"Complete: {len(rows)} null counts")
        print(f"Summary: {output_dir / 'summary.md'}")
        return 0
    except KeyboardInterrupt:
        print("\nInterrupted; completed seeds remain checkpointed.", file=sys.stderr)
        return 130
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
