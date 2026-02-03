#!/usr/bin/env python3
"""
Plot Figure R1: "Stamping becomes the bottleneck (and coloring fixes it)".

Usage:
  python3 color/plot_coloring_benchmark.py
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


SCRIPT_DIR = Path(__file__).parent
DATASET_FILE = SCRIPT_DIR / "coloring_benchmark_results.csv"
NUM_INSTANCES = 1000
NODE_DISTRIBUTION = 39
METHODS = ["loadomp", "loadompColor2"]
TITLE = "Stamping becomes the bottleneck (and coloring fixes it)"


def parse_float(value: str) -> float:
    if value is None:
        return math.nan
    text = value.strip()
    if text == "" or text.lower() == "nan":
        return math.nan
    return float(text)


def parse_int(value: str) -> int:
    if value is None:
        raise ValueError("missing integer value")
    text = value.strip()
    if text == "" or text.lower() == "nan":
        raise ValueError("missing integer value")
    return int(text)


def load_rows(csv_path: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for idx, row in enumerate(reader, start=2):
            try:
                rows.append(
                    {
                        "Method": (row.get("Method") or "").strip(),
                        "NumInstances": parse_int(row.get("NumInstances")),
                        "NumIterations": parse_int(row.get("NumIterations")),
                        "NodeDistributionCount": parse_int(row.get("NodeDistributionCount")),
                        "TargetColors": parse_int(row.get("TargetColors")),
                        "ActualColors": parse_int(row.get("ActualColors")),
                        "ColoringTime_s": parse_float(row.get("ColoringTime_s")),
                        "SingleThreadTime_s": parse_float(row.get("SingleThreadTime_s")),
                        "NumThreads": parse_int(row.get("NumThreads")),
                        "ParallelTime_s": parse_float(row.get("ParallelTime_s")),
                        "ParallelCalcTime_s": parse_float(row.get("ParallelCalcTime_s")),
                        "ParallelApplyTime_s": parse_float(row.get("ParallelApplyTime_s")),
                        "ApplyStamps_pct": parse_float(row.get("ApplyStamps_pct")),
                        "Speedup": parse_float(row.get("Speedup")),
                        "Efficiency_pct": parse_float(row.get("Efficiency_pct")),
                    }
                )
            except ValueError as exc:
                raise ValueError(f"{csv_path}: parse error on line {idx}: {exc}") from exc
    return rows


def is_nan(value: float) -> bool:
    return isinstance(value, float) and math.isnan(value)


def filter_rows(rows: List[Dict[str, object]], method: str) -> List[Dict[str, object]]:
    filtered = [
        r
        for r in rows
        if r["Method"] == method
        and r["NumInstances"] == NUM_INSTANCES
        and r["NodeDistributionCount"] == NODE_DISTRIBUTION
        and r["NumIterations"] > 0
        and not is_nan(r["ParallelCalcTime_s"])
        and not is_nan(r["ParallelApplyTime_s"])
    ]
    if not filtered:
        raise ValueError(
            "No rows found for selection: "
            f"method={method}, NumInstances={NUM_INSTANCES}, "
            f"NodeDistributionCount={NODE_DISTRIBUTION}."
        )
    return filtered


def plot_figure_r1(rows_by_method: Dict[str, List[Dict[str, object]]], out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for ax, method in zip(axes, METHODS):
        rows_sorted = sorted(rows_by_method[method], key=lambda r: r["NumThreads"])
        threads = [r["NumThreads"] for r in rows_sorted]
        calc_ms = [
            (r["ParallelCalcTime_s"] / r["NumIterations"]) * 1000.0 for r in rows_sorted
        ]
        stamp_ms = [
            (r["ParallelApplyTime_s"] / r["NumIterations"]) * 1000.0 for r in rows_sorted
        ]

        x = list(range(len(threads)))
        bar_calc = ax.bar(x, calc_ms, label="Calc time", color="#4C78A8")
        bar_stamp = ax.bar(x, stamp_ms, bottom=calc_ms, label="Stamp time", color="#F58518")

        ax.set_title(method)
        ax.set_xlabel("Threads")
        ax.set_xticks(x)
        ax.set_xticklabels([str(t) for t in threads])
        ax.set_ylim(bottom=0)
        ax.legend()

    axes[0].set_ylabel("Time per iteration (ms)")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Figure R1 (stacked calc/stamp per iteration)."
    )
    parser.add_argument(
        "--out-dir",
        default=str(SCRIPT_DIR / "plots"),
        help="Output directory for plots.",
    )
    parser.add_argument(
        "--format",
        choices=["png", "pdf"],
        default="pdf",
        help="Image format.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = DATASET_FILE
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    rows = load_rows(csv_path)
    rows_by_method = {method: filter_rows(rows, method) for method in METHODS}

    out_name = f"figure_r1_stamping_bottleneck.{args.format}"
    out_path = out_dir / out_name
    plot_figure_r1(rows_by_method, out_path)

    print(
        "Figure R1: "
        f"NumInstances={NUM_INSTANCES}, NodeDistributionCount={NODE_DISTRIBUTION} -> "
        f"{out_path}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
