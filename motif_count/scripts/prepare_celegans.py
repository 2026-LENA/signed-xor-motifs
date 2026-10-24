#!/usr/bin/env python3
"""Convert the Cook et al. SI5 matrix to the C counter's edge-list schema.

The official Cook workbook contains connectivity and weights, but no functional
EXC/INH column.  This converter therefore creates a clearly labelled
neurotransmitter-based approximation using the neuron annotations distributed
with Fenyves et al. (2020). It must not be confused with the separate per-edge
curation used for the paper's reported motif counts.
"""

from __future__ import annotations

import argparse
import csv
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


XML_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def column_number(cell_reference: str) -> int:
    """Return a one-based column number for an Excel cell reference."""
    match = re.match(r"[A-Z]+", cell_reference)
    if match is None:
        raise ValueError(f"Invalid Excel cell reference: {cell_reference}")
    result = 0
    for character in match.group():
        result = result * 26 + ord(character) - ord("A") + 1
    return result


def read_sheet_rows(workbook: Path) -> dict[int, dict[int, str]]:
    """Read the first XLSX sheet with only Python's standard library."""
    with zipfile.ZipFile(workbook) as archive:
        shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        shared_strings = [
            "".join(text.text or "" for text in item.iter(XML_NS + "t"))
            for item in shared_root.findall(XML_NS + "si")
        ]

        sheet_root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        rows: dict[int, dict[int, str]] = {}
        for row in sheet_root.iter(XML_NS + "row"):
            values: dict[int, str] = {}
            for cell in row.findall(XML_NS + "c"):
                value = cell.find(XML_NS + "v")
                if value is None:
                    continue
                text = value.text or ""
                if cell.attrib.get("t") == "s":
                    text = shared_strings[int(text)]
                values[column_number(cell.attrib["r"])] = text
            rows[int(row.attrib["r"])] = values
        return rows


def read_transmitters(path: Path) -> dict[str, str]:
    transmitters: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter=";"):
            transmitter = row["PrimNT"].strip().upper()
            transmitters.setdefault(row["Source"].strip(), transmitter)
    return transmitters


def classification(transmitter: str) -> str:
    if transmitter in {"ACH", "GLU"}:
        return "EXC"
    if transmitter == "GABA":
        return "INH"
    return "OTHER"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("fenyves_csv", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    rows = read_sheet_rows(args.workbook)
    sources = [rows[row_number][3] for row_number in range(4, 304)]
    source_names = set(sources)
    targets = {
        column: name
        for column, name in rows[3].items()
        if 4 <= column <= 456 and name in source_names
    }
    if len(sources) != 300 or len(targets) != 300:
        raise RuntimeError(
            f"Unexpected Cook matrix shape: {len(sources)} sources, "
            f"{len(targets)} targets"
        )

    transmitters = read_transmitters(args.fenyves_csv)
    counts = {"EXC": 0, "INH": 0, "OTHER": 0}
    loop_counts = {"EXC": 0, "INH": 0, "OTHER": 0}
    self_loops = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(
            [
                "Neuron1",
                "Neuron2",
                "Weight",
                "Neurotransmitter",
                "Classification",
                "REF",
                "E/I original",
                "Check",
            ]
        )

        for row_number, source in zip(range(4, 304), sources):
            transmitter = transmitters.get(source, "0")
            edge_class = classification(transmitter)
            for column, target in targets.items():
                raw_weight = rows[row_number].get(column)
                if raw_weight is None or float(raw_weight) == 0.0:
                    continue
                counts[edge_class] += 1
                if source == target:
                    self_loops += 1
                    loop_counts[edge_class] += 1
                writer.writerow(
                    [
                        source,
                        target,
                        int(float(raw_weight)),
                        transmitter,
                        edge_class,
                        "Cook2019-SI5",
                        edge_class,
                        "APPROX_NEURON_NT",
                    ]
                )

    total = sum(counts.values())
    print(f"Wrote {total:,} directed connections to {args.output}")
    print(f"Self-loops (skipped by the C counter): {self_loops:,}")
    print(
        "After self-loop removal: "
        + ", ".join(
            f"{name}={counts[name] - loop_counts[name]:,}"
            for name in ("EXC", "INH", "OTHER")
        )
    )


if __name__ == "__main__":
    main()
