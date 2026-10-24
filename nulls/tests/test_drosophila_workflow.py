#!/usr/bin/env python3
"""Regression tests for Drosophila normalization and null experiments."""

from __future__ import annotations

from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[2]
NORMALIZER = ROOT / "nulls" / "normalize_drosophila.py"
RUNNER = ROOT / "nulls" / "run_drosophila_nulls.py"
RANDOMIZER_SOURCE = ROOT / "nulls" / "randomize_drosophila.c"

HEADER = (
    "pre_root_id",
    "post_root_id",
    "neuropil",
    "syn_count",
    "nt_type",
)


def write_csv(path: Path, rows: list[tuple[object, ...]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(HEADER)
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def run(*command: object, success: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [str(item) for item in command],
        text=True,
        capture_output=True,
        timeout=60,
    )
    if success and result.returncode:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(map(str, command))}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


class DrosophilaWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.build = tempfile.TemporaryDirectory(prefix="drosophila-null-build-")
        cls.randomizer = Path(cls.build.name) / "randomize_drosophila"
        compiler = shutil.which("cc")
        if compiler is None:
            raise unittest.SkipTest("a C compiler is required")
        run(
            compiler,
            "-O2",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            RANDOMIZER_SOURCE,
            "-o",
            cls.randomizer,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.build.cleanup()

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="drosophila-null-test-")
        self.directory = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_normalizer_collapses_signed_edges_and_preserves_synapses(self) -> None:
        source = self.directory / "source.csv"
        output = self.directory / "normalized.csv"
        report = self.directory / "normalization.json"
        write_csv(
            source,
            [
                (1, 2, "A", 3, "ach"),
                (1, 2, "B", 4, "ACH"),
                (1, 2, "C", 5, "GABA"),
                (2, 3, "A", 6, "SER"),
            ],
        )

        run(sys.executable, NORMALIZER, source, output, "--report", report)

        self.assertEqual(
            read_csv(output),
            [
                dict(zip(HEADER, ("1", "2", "AGGREGATED", "7", "ACH"))),
                dict(zip(HEADER, ("1", "2", "AGGREGATED", "5", "GABA"))),
                dict(zip(HEADER, ("2", "3", "AGGREGATED", "6", "SER"))),
            ],
        )
        details = json.loads(report.read_text(encoding="utf-8"))
        self.assertEqual(details["collapse_key"], "nt")
        self.assertEqual(details["original_rows"], 4)
        self.assertEqual(details["unique_edges"], 3)
        self.assertEqual(details["unique_directed_pairs"], 2)
        self.assertEqual(details["collapsed_duplicate_rows"], 1)
        self.assertEqual(details["directed_pairs_with_multiple_output_types"], 1)
        self.assertEqual(details["maximum_rows_collapsed_into_one_edge"], 2)
        self.assertEqual(details["total_syn_count"], 18)

    def test_normalizer_can_collapse_to_binary_signs(self) -> None:
        source = self.directory / "source.csv"
        output = self.directory / "normalized.csv"
        report = self.directory / "normalization.json"
        write_csv(
            source,
            [
                (1, 2, "A", 3, "GABA"),
                (1, 2, "B", 4, "GLUT"),
                (1, 2, "C", 5, "ACH"),
                (2, 3, "A", 6, "SER"),
                (2, 3, "B", 7, "DA"),
            ],
        )

        run(
            sys.executable,
            NORMALIZER,
            source,
            output,
            "--report",
            report,
            "--collapse",
            "sign",
        )

        self.assertEqual(
            read_csv(output),
            [
                dict(zip(HEADER, ("1", "2", "AGGREGATED", "7", "GABA"))),
                dict(zip(HEADER, ("1", "2", "AGGREGATED", "5", "ACH"))),
                dict(zip(HEADER, ("2", "3", "AGGREGATED", "13", "OTHER"))),
            ],
        )
        details = json.loads(report.read_text(encoding="utf-8"))
        self.assertEqual(details["collapse_key"], "sign")
        self.assertEqual(details["unique_edges"], 3)
        self.assertEqual(details["collapsed_duplicate_rows"], 2)
        self.assertEqual(details["total_syn_count"], 25)

    def test_sign_normalizer_rejects_empty_neurotransmitter(self) -> None:
        source = self.directory / "source.csv"
        output = self.directory / "normalized.csv"
        report = self.directory / "normalization.json"
        write_csv(source, [(1, 2, "A", 3, "")])

        result = run(
            sys.executable,
            NORMALIZER,
            source,
            output,
            "--report",
            report,
            "--collapse",
            "sign",
            success=False,
        )

        self.assertEqual(result.returncode, 1)
        self.assertIn("empty pre/post/nt_type", result.stderr)
        self.assertFalse(output.exists())
        self.assertFalse(report.exists())

    def test_normalizer_rejects_colliding_artifact_paths(self) -> None:
        source = self.directory / "source.csv"
        output = self.directory / "normalized.csv"
        write_csv(source, [(1, 2, "A", 3, "ACH")])
        result = run(
            sys.executable,
            NORMALIZER,
            source,
            output,
            "--report",
            output,
            success=False,
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("paths must be different", result.stderr)
        self.assertFalse(output.exists())

    def test_randomizer_preserves_per_type_degrees(self) -> None:
        source = self.directory / "source.csv"
        output = self.directory / "randomized.csv"
        rows = [
            ("a", "b", "A", 1, "ACH"),
            ("c", "d", "A", 1, "ACH"),
            ("e", "f", "A", 1, "ACH"),
            ("g", "h", "A", 1, "ACH"),
            ("i", "j", "A", 1, "GABA"),
            ("k", "l", "A", 1, "GABA"),
            ("m", "n", "A", 1, "GABA"),
            ("o", "p", "A", 1, "GABA"),
            ("q", "r", "A", 1, "SER"),
            ("s", "t", "A", 1, "SER"),
        ]
        write_csv(source, rows)

        run(
            self.randomizer,
            source,
            "--out",
            output,
            "--seed",
            7,
            "--mixing",
            50,
            "--by-nt",
            "-q",
        )

        before = read_csv(source)
        after = read_csv(output)
        self.assertEqual(len(after), len(before))
        self.assertEqual(
            Counter((row["pre_root_id"], row["nt_type"]) for row in after),
            Counter((row["pre_root_id"], row["nt_type"]) for row in before),
        )
        self.assertEqual(
            Counter((row["post_root_id"], row["nt_type"]) for row in after),
            Counter((row["post_root_id"], row["nt_type"]) for row in before),
        )
        self.assertEqual(
            [
                (row["pre_root_id"], row["post_root_id"])
                for row in after
                if row["nt_type"] == "SER"
            ],
            [("q", "r"), ("s", "t")],
        )
        self.assertNotEqual(
            [
                (row["pre_root_id"], row["post_root_id"])
                for row in after
                if row["nt_type"] != "SER"
            ],
            [
                (row["pre_root_id"], row["post_root_id"])
                for row in before
                if row["nt_type"] != "SER"
            ],
        )

    def test_lumped_randomizer_preserves_signed_degrees(self) -> None:
        source = self.directory / "source.csv"
        output = self.directory / "randomized.csv"
        rows = [
            ("a", "b", "A", 1, "ACH"),
            ("c", "d", "A", 1, "ACH"),
            ("e", "f", "A", 1, "ACH"),
            ("g", "h", "A", 1, "ACH"),
            ("i", "j", "A", 1, "GABA"),
            ("k", "l", "A", 1, "GABA"),
            ("m", "n", "A", 1, "GLUT"),
            ("o", "p", "A", 1, "GLUT"),
            ("q", "r", "A", 1, "SER"),
            ("s", "t", "A", 1, "SER"),
        ]
        write_csv(source, rows)

        run(
            self.randomizer,
            source,
            "--out",
            output,
            "--seed",
            11,
            "--mixing",
            50,
            "--lump-inhibitory",
            "-q",
        )

        def sign(row: dict[str, str]) -> str:
            if row["nt_type"] == "ACH":
                return "EXC"
            if row["nt_type"] in {"GABA", "GLUT"}:
                return "INH"
            return "OTHER"

        before = read_csv(source)
        after = read_csv(output)
        self.assertEqual(
            Counter((row["pre_root_id"], sign(row)) for row in after),
            Counter((row["pre_root_id"], sign(row)) for row in before),
        )
        self.assertEqual(
            Counter((row["post_root_id"], sign(row)) for row in after),
            Counter((row["post_root_id"], sign(row)) for row in before),
        )
        self.assertEqual(
            Counter((row["pre_root_id"], row["nt_type"]) for row in after),
            Counter((row["pre_root_id"], row["nt_type"]) for row in before),
        )
        signed_edges = [
            (row["pre_root_id"], row["post_root_id"], sign(row))
            for row in after
            if sign(row) != "OTHER"
        ]
        self.assertEqual(len(signed_edges), len(set(signed_edges)))
        self.assertEqual(
            [
                (row["pre_root_id"], row["post_root_id"])
                for row in after
                if sign(row) == "OTHER"
            ],
            [("q", "r"), ("s", "t")],
        )
        self.assertNotEqual(
            [
                (row["pre_root_id"], row["post_root_id"])
                for row in after
                if sign(row) != "OTHER"
            ],
            [
                (row["pre_root_id"], row["post_root_id"])
                for row in before
                if sign(row) != "OTHER"
            ],
        )

    def test_runner_checkpoints_and_extends_an_experiment(self) -> None:
        source = self.directory / "normalized.csv"
        write_csv(source, [(1, 2, "AGGREGATED", 1, "ACH")])
        fake_randomizer = self.directory / "fake_randomizer.py"
        fake_counter = self.directory / "fake_counter.py"
        fake_randomizer.write_text(
            "#!/usr/bin/env python3\n"
            "from pathlib import Path\n"
            "import shutil, sys\n"
            "args = sys.argv[1:]\n"
            "shutil.copyfile(Path(args[0]), Path(args[args.index('--out') + 1]))\n"
            "print('randomized')\n",
            encoding="utf-8",
        )
        fake_counter.write_text(
            "#!/usr/bin/env python3\n"
            "from pathlib import Path\n"
            "import re, sys\n"
            "match = re.search(r'null_(\\d+)\\.csv$', Path(sys.argv[-1]).name)\n"
            "count = 90 + int(match.group(1)) if match else 100\n"
            "print(f'MOTIF_COUNT_RAW {count}  MOTIF_COUNT_DISTINCT {count}')\n",
            encoding="utf-8",
        )
        fake_randomizer.chmod(0o755)
        fake_counter.chmod(0o755)
        results = self.directory / "results"
        base = (
            sys.executable,
            RUNNER,
            source,
            "--mixing",
            1,
            "--randomizer",
            fake_randomizer,
            "--counter",
            fake_counter,
            "--output-dir",
            results,
        )

        first = run(*base, "--nulls", 2)
        self.assertIn("Complete: 2 null counts", first.stdout)
        first_log = results / "logs" / "null_1_count.log"
        first_mtime = first_log.stat().st_mtime_ns
        time.sleep(0.01)

        second = run(*base, "--nulls", 3)
        self.assertIn("seed 1: already complete", second.stdout)
        self.assertIn("seed 2: already complete", second.stdout)
        self.assertIn("Complete: 3 null counts", second.stdout)
        self.assertEqual(first_log.stat().st_mtime_ns, first_mtime)
        with (results / "null_counts.csv").open(
            encoding="utf-8", newline=""
        ) as stream:
            counts = list(csv.DictReader(stream))
        self.assertEqual(
            [(row["seed"], row["count"]) for row in counts],
            [("1", "91"), ("2", "92"), ("3", "93")],
        )
        summary = json.loads(
            (results / "summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(summary["observed"], 100)
        self.assertEqual(summary["n_nulls"], 3)
        self.assertEqual(summary["null_mean"], 92.0)
        self.assertFalse(any((results / "work").iterdir()))


if __name__ == "__main__":
    unittest.main()
