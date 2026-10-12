#!/usr/bin/env python3
"""CPU reference checks, plus optional comparisons with the selected backend.

Fixtures are synthetic; no downloaded dataset or third-party Python package is
needed. CSV comparisons use multisets because OpenMP output order is arbitrary.
"""
import argparse
from collections import Counter
import csv
from pathlib import Path
import random
import re
import subprocess
import tempfile
import unittest

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--bin-dir", required=True, type=Path)
parser.add_argument("--gpu", action="store_true")
args, unittest_args = parser.parse_known_args()
BIN = args.bin_dir.resolve()
PROGRAMS = ("drosophila_signed_xor_motif", "celegans_signed_xor_motif")


def variant_edges(variant="B", sign5=1, offset=0):
    # n1,n2,n3,n4,n5,n6,F+,F- = 0..7; 8 and 9 label the feedback nodes
    edges = [(0, 1, 1), (0, 5, 1), (2, 3, 1), (2, 5, 1),
             (1, 4, 1), (3, 4, 1), (5, 1, 2), (5, 3, 2),
             (4, 6, sign5), (4, 7, sign5), (6, 8, 1), (7, 9, 2)]
    if "A" in variant:
        edges += [(0, 6, 1), (2, 7, 1)]
    if "B" in variant:
        edges += [(1, 6, 1), (3, 7, 1)]
    return [(s + offset, t + offset, k) for s, t, k in edges]


def write_edges(path, program, edges, labels=None, weights=None):
    with path.open("w", newline="") as fp:
        writer = csv.writer(fp, delimiter=";" if program.startswith("celegans") else ",")
        if program == "signed_xor_motif_count":
            writer.writerow(["source_node_id", "target_node_id", "source_label", "connection_type", "syn_weight"])
            for i, (s, t, k) in enumerate(edges):
                label = (labels or {}).get(s, "i23Pvalb" if s == 5 else "i4Sst" if s == 7 else "e23")
                writer.writerow([s, t, label, "EXC" if k == 1 else "INH", (weights or {}).get(i, 1 if k == 1 else -1)])
        elif program.startswith("celegans"):
            writer.writerow(["Neuron1", "Neuron2", "Weight", "Neurotransmitter", "Classification"])
            for s, t, k in edges:
                writer.writerow([s, t, 1, "ACH" if k == 1 else "GABA", "EXC" if k == 1 else "INH" if k == 2 else "OTHER"])
        else:
            writer.writerow(["pre_root_id", "post_root_id", "neuropil", "syn_count", "nt_type"])
            for s, t, k in edges:
                writer.writerow([s, t, "test", 1, "ACH" if k == 1 else "GABA" if k == 2 else "DA"])


def run(program, path, backend, *options, success=True):
    result = subprocess.run([str(BIN / program), str(path), f"--{backend}", *map(str, options)],
                            text=True, capture_output=True, timeout=120)
    if success and result.returncode:
        raise AssertionError(f"{program} --{backend} failed:\n{result.stderr}\n{result.stdout}")
    return result


def count(program, output):
    pattern = r"^TOTAL\s+(\d+)" if program == "signed_xor_motif_count" else r"MOTIF_COUNT_RAW (\d+)"
    return int(re.search(pattern, output, re.MULTILINE).group(1))


def csv_rows(path):
    with path.open(newline="") as fp:
        return list(csv.DictReader(fp))


class MotifTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="motif-gpu-test-")
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def compare(self, program, edges, expected=None, induced=False, export=True,
                labels=None, weights=None, extra=()):
        path = self.root / "edges.csv"
        write_edges(path, program, edges, labels, weights)
        options = ["--induced"] if induced else []
        options += list(extra)
        if program != "signed_xor_motif_count":
            options += ["-q"]
        cpu_csv = self.root / "cpu.csv"
        if export and program != "signed_xor_motif_count":
            options += ["--out", cpu_csv]
        cpu = run(program, path, "cpu", *options)
        cpu_count = count(program, cpu.stdout)
        if expected is not None:
            self.assertEqual(cpu_count, expected, cpu.stdout)
        if export and program != "signed_xor_motif_count":
            self.assertEqual(len(csv_rows(cpu_csv)), cpu_count)
        if args.gpu:
            gpu_csv = self.root / "gpu.csv"
            gpu_options = [gpu_csv if x == cpu_csv else x for x in options]
            gpu = run(program, path, "gpu", *gpu_options)
            self.assertEqual(count(program, gpu.stdout), cpu_count)
            if export and program != "signed_xor_motif_count":
                cpu_rows = Counter(tuple(sorted(row.items())) for row in csv_rows(cpu_csv))
                gpu_rows = Counter(tuple(sorted(row.items())) for row in csv_rows(gpu_csv))
                self.assertEqual(gpu_rows, cpu_rows)
            if program == "signed_xor_motif_count":
                # Check each inhibitory type, not just the grand total.
                self.assertEqual(re.findall(r"^i\S+\s+\d+\s+\d+", gpu.stdout, re.MULTILINE),
                                 re.findall(r"^i\S+\s+\d+\s+\d+", cpu.stdout, re.MULTILINE))
        return cpu_count

    def test_feedback_variants_and_xor_sign(self):
        for program in PROGRAMS:
            for variant in ("A", "B", "AB"):
                for sign in (1, 2):
                    for induced in (False, True):
                        with self.subTest(program=program, variant=variant, sign=sign, induced=induced):
                            self.compare(program, variant_edges(variant, sign), 1, induced)

    def test_induced_extra_wrong_sign_and_other_edges(self):
        for program in PROGRAMS:
            for edge in ((6, 0, 1), (0, 1, 2), (0, 1, 3)):
                with self.subTest(program=program, edge=edge):
                    self.compare(program, variant_edges() + [edge], induced=False)
                    self.compare(program, variant_edges() + [edge], 0, induced=True)

    def test_feedback_identity(self):
        for program in PROGRAMS:
            for edges in (variant_edges() + [(6, 9, 2)],
                          [e for e in variant_edges() if e[0] != 6],
                          variant_edges() + [(4, 9, 2)],
                          variant_edges() + [(7, 8, 1)]):
                self.compare(program, edges, 0)

    def test_feedback_inclusion_exclusion(self):
        edges = variant_edges("AB") + [(4, 10, 1), (4, 11, 1), (4, 12, 1), (4, 13, 1),
                                        (10, 8, 1), (11, 8, 1), (12, 9, 2), (13, 9, 2),
                                        (0, 10, 1), (1, 11, 1), (2, 12, 1), (3, 13, 1)]
        for program in PROGRAMS:
            self.compare(program, edges, 7)
            self.compare(program, edges, 3, induced=True)

    def test_empty_missing_core_and_distinctness(self):
        for program in (*PROGRAMS, "signed_xor_motif_count"):
            for edges in ([], [(0, 1, 1)], variant_edges()[1:],
                          [(0 if s == 7 else s, 0 if t == 7 else t, k)
                           for s, t, k in variant_edges()]):
                self.compare(program, edges, 0)

    def test_random_graphs(self):
        rng = random.Random(55123)
        for seed in range(5):
            edges = variant_edges("AB", 1 if seed % 2 else 2)
            for _ in range(15):
                s, t = rng.randrange(12), rng.randrange(12)
                edges.append((s, t, rng.choice((1, 2, 3))))
            rng.shuffle(edges)
            for program in PROGRAMS:
                self.compare(program, edges)
                self.compare(program, edges, induced=True)
            # Allen input has signed edges only.
            self.compare("signed_xor_motif_count", [e for e in edges if e[2] != 3])

    def test_duplicate_edges_and_output_pagination(self):
        edges = variant_edges()
        edges += [(4, 6, 1)] * 299 + [(4, 7, 1)] * 299
        for program in (*PROGRAMS, "signed_xor_motif_count"):
            self.compare(program, edges, 90000)

    def test_multiple_task_batches_and_input_slices(self):
        edges = variant_edges() + [(0, 1, 1)] * 4199 + [(2, 3, 1)] * 8
        for program in (*PROGRAMS, "signed_xor_motif_count"):
            self.compare(program, edges, 37800)

    def test_batch_boundary(self):
        # Each duplicate n1->n2 yields one job. Cross the 32,768-job batch
        # boundary without creating a large motif output file.
        edges = variant_edges() + [(0, 1, 1)] * 32770
        for program in PROGRAMS:
            self.compare(program, edges, 32771, export=False)

    def test_allen_order_limits_pivot_and_weights(self):
        program = "signed_xor_motif_count"
        edges = variant_edges()
        self.compare(program, edges, 1)
        # The CPU uses i<j in the original INH adjacency order. Sorting that
        # iteration order would change which target receives F+ versus F-.
        swapped = [e for e in edges if e[0] != 5] + [(5, 3, 2), (5, 1, 2)]
        self.compare(program, swapped, 0)
        self.compare(program, edges, 1, extra=("--pivot", 5, "-n", 1))
        self.compare(program, edges, 0, weights={0: 0.01}, extra=("--min-weight", 0.1))
        self.compare(program, edges, 1, weights={6: -0.5}, extra=("--min-weight", 0.1))
        # Give F- another target type to check per-type reporting.
        self.compare(program, edges, 1, labels={5: "i6Sst", 7: "i23Pvalb"})
        path = self.root / "edges.csv"
        weights = {i: (i + 0.5) * (1 if k == 1 else -1) for i, (_, _, k) in enumerate(edges)}
        write_edges(path, program, edges, weights=weights)
        for backend in ("cpu", "gpu") if args.gpu else ("cpu",):
            prefix = self.root / backend
            output = run(program, path, backend, "--weight-stats", prefix, "--sample-size", 10)
            self.assertEqual(count(program, output.stdout), 1)
            stats = csv_rows(Path(str(prefix) + "_stats.csv"))
            self.assertEqual(len(stats), 12)
            self.assertTrue(all(row["n_observations"] == "1" for row in stats))
            self.assertEqual(len(csv_rows(Path(str(prefix) + "_sample.csv"))), 12)
        if args.gpu:
            for suffix in ("_stats.csv", "_histogram.csv", "_sample.csv"):
                self.assertEqual(Path(str(self.root / "cpu") + suffix).read_text(),
                                 Path(str(self.root / "gpu") + suffix).read_text())

    def test_paginated_allen_weight_stats(self):
        program = "signed_xor_motif_count"
        edges = variant_edges() + [(4, 6, 1)] * 299 + [(4, 7, 1)] * 299
        path = self.root / "edges.csv"
        write_edges(path, program, edges)
        for backend in ("cpu", "gpu") if args.gpu else ("cpu",):
            prefix = self.root / ("many-" + backend)
            result = run(program, path, backend, "--weight-stats", prefix, "--sample-size", 10)
            self.assertEqual(count(program, result.stdout), 90000)
            stats = csv_rows(Path(str(prefix) + "_stats.csv"))
            self.assertTrue(all(row["n_observations"] == "90000" for row in stats))
            self.assertEqual(len(csv_rows(Path(str(prefix) + "_sample.csv"))), 120)
        if args.gpu:
            for suffix in ("_stats.csv", "_histogram.csv"):
                self.assertEqual(Path(str(self.root / "many-cpu") + suffix).read_text(),
                                 Path(str(self.root / "many-gpu") + suffix).read_text())

    def test_gpu_flag_on_cpu_only_build(self):
        if args.gpu:
            return
        path = self.root / "edges.csv"
        for program in (*PROGRAMS, "signed_xor_motif_count"):
            write_edges(path, program, variant_edges())
            result = run(program, path, "gpu", success=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("GPU support is not compiled in", result.stderr)


if __name__ == "__main__":
    unittest.main(argv=[__file__, *unittest_args])
