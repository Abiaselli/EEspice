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
DATASET_FILE_10K = SCRIPT_DIR / "coloring_benchmark_results10k.csv"
NUM_INSTANCES = 1000
NODE_DISTRIBUTION = 39
METHODS = ["loadomp", "loadompColor"]
TITLE = "Stamping becomes the bottleneck (and coloring fixes it)"

# Figure R2 configuration
R2_NUM_INSTANCES = 1000
R2_NODE_DISTRIBUTIONS = [
    {"node_dist": 1000, "label": "(a) Low conflict"},  # Keep hard-coded for low conflict
    {"node_dist_selector": "moderate", "label": "(b) Moderate conflict"},
    {"node_dist_selector": "min", "label": "(c) High conflict"},
]
R2_METHODS = ["loadomp", "loadompColor", "loadompColorFused"]

# Figure R3 configuration
R3_NUM_THREADS = 64  # Fixed thread count
R3_METHODS = ["loadomp", "loadompColor", "loadompColorFused"]
R3_SUBPLOTS = [
    {"num_instances": 100, "label": "(a) 100 instances"},
    {"num_instances": 1000, "label": "(b) 1000 instances"},
    {"num_instances": 10000, "label": "(c) 10000 instances"},
]

# Figure R4 configuration
R4_METHODS = ["loadomp", "loadompColor", "loadompColorFused"]
R4_CONFLICT_LEVELS = [
    {"level": "low", "label": "Low conflict", "node_dist_selector": "equal"},
    {"level": "moderate", "label": "Moderate conflict", "node_dist_selector": "moderate"},
    {"level": "high", "label": "High conflict", "node_dist_selector": "min"},  # smallest available
]
R4_NUM_INSTANCES = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                    2000, 4000, 6000, 8000, 10000]

METHOD_STYLES = {
    "loadomp": {"color": "#E45756", "marker": "s", "label": "loadomp"},
    "loadompColor": {"color": "#4C78A8", "marker": "o", "label": "loadompColor"},
    "loadompColorFused": {"color": "#72B7B2", "marker": "^", "label": "loadompColorFused"},
}


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


def filter_rows_r2(
    rows: List[Dict[str, object]], method: str, node_distribution: int
) -> List[Dict[str, object]]:
    """Filter rows for Figure R2 by method and NodeDistributionCount."""
    return [
        r
        for r in rows
        if r["Method"] == method
        and r["NumInstances"] == R2_NUM_INSTANCES
        and r["NodeDistributionCount"] == node_distribution
        and not is_nan(r["Speedup"])
    ]


def filter_rows_r3(
    rows: List[Dict[str, object]], method: str, num_threads: int, num_instances: int
) -> List[Dict[str, object]]:
    """Filter rows for Figure R3: fixed instance count and thread count."""
    return [
        r
        for r in rows
        if r["Method"] == method
        and r["NumInstances"] == num_instances
        and r["NumThreads"] == num_threads
        and not is_nan(r["Speedup"])
    ]


def get_available_node_dists(
    rows: List[Dict[str, object]], num_instances: int, method: str = None
) -> List[int]:
    """Get all unique NodeDistributionCount values for a given NumInstances and optional method."""
    return sorted(set(
        r["NodeDistributionCount"]
        for r in rows
        if r["NumInstances"] == num_instances
        and (method is None or r["Method"] == method)
    ))


def get_node_dist_for_conflict_level(
    rows: List[Dict[str, object]], num_instances: int, selector: object, method: str = None
) -> int:
    """
    Select appropriate NodeDistributionCount based on conflict level.
    - "equal": Returns NumInstances (low conflict, ~1 color)
    - "moderate": Returns closest available value around 40 (~N/40 effective colors)
    - "min": Returns smallest available NodeDistributionCount (~N/3 colors)
    - int: Returns that specific value (e.g., 3 for high conflict)
    """
    available = get_available_node_dists(rows, num_instances, method)
    if not available:
        return -1

    if selector == "equal":
        # Low conflict: NodeDistributionCount = NumInstances
        if num_instances in available:
            return num_instances
        return max(available)  # fallback to largest
    elif selector == "moderate":
        # Moderate conflict: target around 40 (~N/40 effective colors)
        target = 40
        closest = min(available, key=lambda x: abs(x - target))
        return closest
    elif selector == "min":
        # High conflict: smallest available NodeDistributionCount
        return min(available)
    elif isinstance(selector, int):
        # High conflict: specific value (e.g., 3)
        if selector in available:
            return selector
        # Fallback to closest available
        closest = min(available, key=lambda x: abs(x - selector))
        return closest
    return available[0]


def filter_rows_r4(
    rows: List[Dict[str, object]], method: str, num_instances: int, node_distribution: int
) -> List[Dict[str, object]]:
    """Filter rows for Figure R4 by method, num_instances, and node_distribution."""
    return [
        r
        for r in rows
        if r["Method"] == method
        and r["NumInstances"] == num_instances
        and r["NodeDistributionCount"] == node_distribution
        and not is_nan(r["Speedup"])
    ]


def compute_effective_colors(row: Dict[str, object]) -> int:
    """Compute effective colors: ActualColors if > 0, else ceil(NumInstances/NodeDistributionCount)."""
    if row["ActualColors"] > 0:
        return row["ActualColors"]
    return math.ceil(row["NumInstances"] / row["NodeDistributionCount"])


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


def plot_figure_r2(rows: List[Dict[str, object]], out_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    for ax, config in zip(axes, R2_NODE_DISTRIBUTIONS):
        label = config["label"]

        # Resolve node_dist: use hard-coded value or dynamic selector
        if "node_dist" in config:
            node_dist = config["node_dist"]
        else:
            selector = config["node_dist_selector"]
            node_dist = get_node_dist_for_conflict_level(rows, R2_NUM_INSTANCES, selector)

        # Track ActualColors to compute subtitle dynamically
        actual_colors_set = set()

        for method in R2_METHODS:
            style = METHOD_STYLES[method]
            filtered = filter_rows_r2(rows, method, node_dist)
            if not filtered:
                continue
            rows_sorted = sorted(filtered, key=lambda r: r["NumThreads"])
            threads = [r["NumThreads"] for r in rows_sorted]
            speedups = [r["Speedup"] for r in rows_sorted]

            # Collect ActualColors from filtered data
            for r in filtered:
                actual_colors_set.add(r["ActualColors"])

            ax.plot(
                threads,
                speedups,
                marker=style["marker"],
                color=style["color"],
                label=style["label"],
                linewidth=1.5,
                markersize=5,
            )

        # Compute subtitle from actual data
        if actual_colors_set:
            actual_colors_val = max(actual_colors_set)  # Use max if multiple values
            subtitle = f"ActualColors = {actual_colors_val}"
        else:
            subtitle = ""

        ax.set_title(f"{label}\n{subtitle}", fontsize=10)
        ax.set_xlabel("Threads")
        ax.set_ylabel("Speedup")
        ax.set_xscale("log", base=2)
        ax.set_xticks([2, 4, 8, 16, 32, 64, 128])
        ax.set_xticklabels(["2", "4", "8", "16", "32", "64", "128"])
        ax.set_ylim(bottom=0)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02))

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_figure_r3(
    rows: List[Dict[str, object]], rows_10k: List[Dict[str, object]], out_path: Path
) -> None:
    """
    Figure R3: Speedup vs EffectiveColors at fixed thread count.
    Three subplots for different instance counts (100, 1000, 10000).
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    for ax, config in zip(axes, R3_SUBPLOTS):
        num_instances = config["num_instances"]
        label = config["label"]

        # Select appropriate data source
        source_rows = rows_10k if num_instances == 10000 else rows

        for method in R3_METHODS:
            method_rows = filter_rows_r3(source_rows, method, R3_NUM_THREADS, num_instances)
            if not method_rows:
                continue

            # Compute effective colors and sort
            data = [(compute_effective_colors(r), r["Speedup"]) for r in method_rows]
            data.sort(key=lambda x: x[0])

            colors_list = [d[0] for d in data]
            speedups = [d[1] for d in data]

            style = METHOD_STYLES[method]
            ax.plot(
                colors_list,
                speedups,
                marker=style["marker"],
                color=style["color"],
                label=style["label"],
                linewidth=1.5,
                markersize=5,
            )

        ax.set_xlabel("Effective Colors")
        ax.set_ylabel("Speedup")
        ax.set_title(label)
        ax.set_xscale("log")
        ax.set_yscale("log", base=2)

    # Shared legend at top
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02))

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_figure_r4(
    rows: List[Dict[str, object]], rows_10k: List[Dict[str, object]], out_path: Path
) -> None:
    """
    Figure R4: 3×3 grid showing Speedup vs Threads.
    Rows: Methods (loadomp, loadompColor, loadompColorFused)
    Columns: Conflict levels (low, moderate, high)
    Lines: Different NumInstances values with color gradient.
    """
    import matplotlib.cm as cm

    # Each subplot has independent y-axis for best data visualization
    fig, axes = plt.subplots(3, 3, figsize=(16, 12), sharex=True)

    # Combine data from both sources
    all_rows = rows + rows_10k

    # Create color map for NumInstances values with shared norm for consistency
    num_instances_list = R4_NUM_INSTANCES
    cmap = cm.viridis
    norm = plt.Normalize(vmin=min(num_instances_list), vmax=max(num_instances_list))

    for row_idx, method in enumerate(R4_METHODS):
        for col_idx, conflict_config in enumerate(R4_CONFLICT_LEVELS):
            ax = axes[row_idx, col_idx]
            selector = conflict_config["node_dist_selector"]

            for num_instances in num_instances_list:
                # Get appropriate node distribution for this conflict level (method-aware)
                node_dist = get_node_dist_for_conflict_level(
                    all_rows, num_instances, selector, method=method
                )
                if node_dist < 0:
                    continue

                filtered = filter_rows_r4(all_rows, method, num_instances, node_dist)
                if not filtered:
                    continue

                rows_sorted = sorted(filtered, key=lambda r: r["NumThreads"])
                threads = [r["NumThreads"] for r in rows_sorted]
                speedups = [r["Speedup"] for r in rows_sorted]

                ax.plot(
                    threads,
                    speedups,
                    marker="o",
                    color=cmap(norm(num_instances)),
                    linewidth=1.2,
                    markersize=4,
                    alpha=0.8,
                )

            # Set column titles only on top row
            if row_idx == 0:
                ax.set_title(f"{conflict_config['label']}", fontsize=11)

            # Set x-axis label only on bottom row
            if row_idx == 2:
                ax.set_xlabel("Threads")

            # Set y-axis label only on left column
            if col_idx == 0:
                ax.set_ylabel("Speedup")

            ax.set_xscale("log", base=2)
            ax.set_xticks([1, 2, 4, 8, 16, 32, 64, 128])
            ax.set_xticklabels(["1", "2", "4", "8", "16", "32", "64", "128"])

    # Adjust layout for row labels on left and colorbar on right
    fig.subplots_adjust(left=0.14, right=0.88)

    # Add row labels (method names) dynamically based on subplot positions
    for row_idx, method in enumerate(R4_METHODS):
        pos = axes[row_idx, 0].get_position()
        y_center = 0.5 * (pos.y0 + pos.y1)
        x_text = pos.x0 - 0.06
        fig.text(
            x_text,
            y_center,
            method,
            fontsize=11,
            fontweight="bold",
            rotation=90,
            va="center",
            ha="center",
        )

    # Create colorbar for NumInstances on the right side (using same norm as lines)
    cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="vertical")
    cbar.set_label("NumInstances", fontsize=10)

    fig.suptitle("Speedup vs Threads (by Method and Conflict Level)", fontsize=14, y=0.98)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot benchmark figures (R1: stacked calc/stamp, R2: speedup vs threads, R3: speedup vs colors, R4: speedup grid)."
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
    parser.add_argument(
        "--figure",
        choices=["r1", "r2", "r3", "r4", "all"],
        default="all",
        help="Which figure(s) to generate.",
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

    if args.figure in ("r1", "all"):
        rows_by_method = {method: filter_rows(rows, method) for method in METHODS}
        out_name = f"figure_r1_stamping_bottleneck.{args.format}"
        out_path = out_dir / out_name
        plot_figure_r1(rows_by_method, out_path)
        print(
            "Figure R1: "
            f"NumInstances={NUM_INSTANCES}, NodeDistributionCount={NODE_DISTRIBUTION} -> "
            f"{out_path}"
        )

    if args.figure in ("r2", "all"):
        out_name = f"figure_r2_speedup_vs_threads.{args.format}"
        out_path = out_dir / out_name
        plot_figure_r2(rows, out_path)
        node_dist_info = [
            c.get("node_dist") or c.get("node_dist_selector") for c in R2_NODE_DISTRIBUTIONS
        ]
        print(
            "Figure R2: "
            f"NumInstances={R2_NUM_INSTANCES}, "
            f"NodeDistributions={node_dist_info} -> "
            f"{out_path}"
        )

    if args.figure in ("r3", "all"):
        # Load 10k dataset for R3
        csv_path_10k = DATASET_FILE_10K
        if not csv_path_10k.exists():
            raise FileNotFoundError(f"CSV not found: {csv_path_10k}")
        rows_10k = load_rows(csv_path_10k)

        out_name = f"figure_r3_speedup_vs_colors.{args.format}"
        out_path = out_dir / out_name
        plot_figure_r3(rows, rows_10k, out_path)
        print(
            "Figure R3: "
            f"NumInstances={[c['num_instances'] for c in R3_SUBPLOTS]}, "
            f"NumThreads={R3_NUM_THREADS} -> "
            f"{out_path}"
        )

    if args.figure in ("r4", "all"):
        # Load 10k dataset for R4
        csv_path_10k = DATASET_FILE_10K
        if not csv_path_10k.exists():
            raise FileNotFoundError(f"CSV not found: {csv_path_10k}")
        rows_10k = load_rows(csv_path_10k)

        out_name = f"figure_r4_speedup_grid.{args.format}"
        out_path = out_dir / out_name
        plot_figure_r4(rows, rows_10k, out_path)
        print(
            "Figure R4: "
            f"3x3 grid (Methods x Conflict Levels), "
            f"NumInstances={R4_NUM_INSTANCES} -> "
            f"{out_path}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
