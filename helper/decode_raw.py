#!/usr/bin/env python3
"""Minimal decoder for EEspice ngspice-compatible binary .raw files.

The C++ writer emits native binary doubles. EEspice's supported Linux/x86_64
output is little-endian, so this example decodes values with "<d" and "<dd".
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path


def read_raw_plot(f):
    """Read one RAW plot block from an open binary file, or return None at EOF."""
    header = {}
    variables = []

    first = f.readline()
    if not first:
        return None
    while first == b"\n":
        first = f.readline()
        if not first:
            return None

    key, value = first.decode("ascii").rstrip("\n").split(": ", 1)
    header[key] = value

    while True:
        line = f.readline().decode("ascii").rstrip("\n")
        if line == "Variables:":
            break
        key, value = line.split(": ", 1)
        header[key] = value

    for _ in range(int(header["No. Variables"])):
        index, name, typ = f.readline().decode("ascii").strip().split("\t")
        variables.append((int(index), name, typ))

    if f.readline() != b"Binary:\n":
        raise ValueError("RAW plot is missing Binary delimiter")

    nvars = int(header["No. Variables"])
    npoints = int(header["No. Points"])
    is_complex = header["Flags"] == "complex"
    item_size = 16 if is_complex else 8
    blob = f.read(npoints * nvars * item_size)

    rows = []
    off = 0
    for _ in range(npoints):
        row = []
        for _ in range(nvars):
            if is_complex:
                re, im = struct.unpack_from("<dd", blob, off)
                row.append(complex(re, im))
            else:
                (value,) = struct.unpack_from("<d", blob, off)
                row.append(value)
            off += item_size
        rows.append(row)

    return header, variables, rows


def read_raw_file(path):
    """Read every plot block from a RAW file."""
    plots = []
    with open(path, "rb") as f:
        while True:
            plot = read_raw_plot(f)
            if plot is None:
                break
            plots.append(plot)
    return plots


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Decode an EEspice binary .raw file and print plot summaries."
    )
    parser.add_argument("raw_file", type=Path)
    args = parser.parse_args()

    for plot_index, (header, variables, rows) in enumerate(
        read_raw_file(args.raw_file), start=1
    ):
        names = [name for _, name, _ in variables]
        print(
            f"plot {plot_index}: {header['Plotname']} "
            f"{header['Flags']} {len(variables)} vars {len(rows)} points"
        )
        print("variables:", ", ".join(names))
        if rows:
            print("first row:", rows[0])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
