#!/usr/bin/env python3
"""
Cache Miss Analysis Script

This script runs danial_cache.o with varying array sizes and collects
L1, L2, and L3 cache miss statistics using perf. It then generates plots
showing the relationship between array size and cache misses.

L3 miss rate is calculated as cache-misses / cache-references.
"""

import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Configuration
EXECUTABLE = "./danial_cache.o"
ARRAY_SIZES = [
    1000, 2000, 4000, 8000,           # < 64KB (fits in L1)
    16000, 32000, 64000,              # ~128KB-512KB (fits in L2)
    128000, 256000, 512000,           # ~1MB-4MB (exceeds L2, fits in L3)
    1000000, 2000000, 4000000         # > 8MB (exceeds most caches)
]
ITERATION_CONFIGS = [1, 100000]  # Test with 1 iteration and 100000 iterations

# Perf events to monitor
# L1: L1-dcache-loads and L1-dcache-load-misses for accurate L1 miss rate
# L3 (LLC): Generic cache-references and cache-misses (Last Level Cache)
#           L3 miss rate = cache-misses / cache-references
PERF_EVENTS = "L1-dcache-loads,L1-dcache-load-misses,cache-references,cache-misses"

# AMD-specific L2 metrics (must use -M flag, may not be available on all kernels)
# Note: These may not work without elevated permissions
PERF_METRICS = "all_l2_cache_accesses,all_l2_cache_misses"


def run_perf_test(array_size, iterations):
    """
    Run danial_cache.o with perf stat to collect cache statistics.

    Args:
        array_size: Size of array to test
        iterations: Number of iterations to run

    Returns:
        Dictionary with cache statistics or None if failed
    """
    # Build perf command
    events = PERF_EVENTS

    cmd = [
        "perf", "stat",
        "--no-big-num",            # Disable locale/group separators
        "-x", ";",                 # CSV-like output with semicolon delimiter
        "-e", events,
        "-M", PERF_METRICS,
        EXECUTABLE,
        str(array_size),
        str(iterations)
    ]

    try:
        # Run perf stat (output goes to stderr)
        # Set LC_ALL=C to ensure consistent decimal formatting (dot separator)
        env = dict(os.environ, LC_ALL="C")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60*10,  # 10 minutes timeout
            env=env
        )

        # Parse perf output
        stats = parse_perf_output(result.stderr)
        stats['array_size'] = array_size
        stats['iterations'] = iterations
        stats['array_size_kb'] = (array_size * 8) / 1024  # doubles are 8 bytes

        # Diagnostic: warn if all counters are zero (likely parsing or permission issue)
        all_zero = all(stats[key] == 0 for key in ['l1_loads', 'l1_misses', 'l2_accesses', 'l2_misses', 'l3_accesses', 'l3_misses'])
        if all_zero:
            print(f"\n WARNING: All counters are zero for array_size={array_size}, iterations={iterations}")
            print("    This may indicate a parsing or permission issue.")
            print("    Raw perf output:")
            for line in result.stderr.splitlines()[:20]:  # Print first 20 lines
                print(f"    {line}")
            if len(result.stderr.splitlines()) > 20:
                print(f"    ... ({len(result.stderr.splitlines()) - 20} more lines)")
            print()

        return stats

    except subprocess.TimeoutExpired:
        print(f"Timeout for size={array_size}, iterations={iterations}")
        return None
    except Exception as e:
        print(f"Error running perf for size={array_size}: {e}")
        return None


def parse_perf_output(output: str):
    """
    Parse perf stat CSV output to extract cache statistics.

    With --no-big-num and -x ";", perf stat outputs CSV format:
    - Event lines: value;unit;event;runtime;pcnt;[var];[metric_value];[metric_unit]
    - Metric lines: Additional metrics with empty earlier fields, only metric_value;metric_name filled

    Args:
        output: stderr output from perf stat (CSV format)

    Returns:
        Dictionary with parsed statistics
    """
    stats = {
        'l1_loads': 0,
        'l1_misses': 0,
        'l2_accesses': 0,
        'l2_misses': 0,
        'l3_accesses': 0,
        'l3_misses': 0
    }

    def safe_float(value_str):
        """Safely convert string to float, return None if invalid."""
        try:
            return float(value_str)
        except (ValueError, TypeError):
            return None

    # Parse CSV lines
    for line in output.splitlines():
        line = line.strip()

        # Skip empty lines and error messages
        if not line or "<not supported>" in line or "<not counted>" in line:
            continue

        parts = [p.strip() for p in line.split(";")]

        # CSV format: value;unit;event;runtime;pcnt;[optional fields];[metric_value];[metric_unit]
        # Metric lines have empty earlier fields and data in the last two columns

        # 1) Try to parse as metric line (metric_value;metric_name in last two columns)
        if len(parts) >= 2:
            metric_value = safe_float(parts[-2])
            metric_name = parts[-1]

            if metric_value is not None and metric_name:
                if metric_name == "all_l2_cache_accesses":
                    stats['l2_accesses'] = int(metric_value)
                    continue
                elif metric_name == "all_l2_cache_misses":
                    stats['l2_misses'] = int(metric_value)
                    continue

        # 2) Try to parse as event line (value;unit;event;...)
        if len(parts) >= 3:
            value = safe_float(parts[0])
            event_name = parts[2]

            if value is None or not event_name:
                continue

            # Match event names
            if event_name == "L1-dcache-loads":
                stats['l1_loads'] = int(value)
            elif event_name == "L1-dcache-load-misses":
                stats['l1_misses'] = int(value)
            elif event_name == "cache-references":
                # Use generic LLC references for L3 accesses
                stats['l3_accesses'] = int(value)
            elif event_name == "cache-misses":
                # Use generic LLC misses for L3 misses
                stats['l3_misses'] = int(value)

    return stats


def collect_data():
    """
    Collect cache miss data for all array sizes and iteration counts.

    Returns:
        pandas DataFrame with all collected data
    """
    all_data = []
    total_tests = len(ARRAY_SIZES) * len(ITERATION_CONFIGS)
    test_count = 0

    print(f"Running {total_tests} tests...")
    print(f"Array sizes: {len(ARRAY_SIZES)} configurations")
    print(f"Iteration counts: {ITERATION_CONFIGS}")
    print("-" * 60)

    for iterations in ITERATION_CONFIGS:
        print(f"\nTesting with {iterations} iteration(s):")
        for array_size in ARRAY_SIZES:
            test_count += 1
            print(f"  [{test_count}/{total_tests}] Array size: {array_size:>8} elements "
                  f"({(array_size*8/1024):.1f} KB)...", end=" ", flush=True)

            stats = run_perf_test(array_size, iterations)

            if stats:
                all_data.append(stats)
                print(f"L1: {stats['l1_misses']:>10,}, L2: {stats['l2_misses']:>10,}, L3: {stats['l3_misses']:>10,}")
            else:
                print("FAILED")

    print("\n" + "=" * 60)
    print(f"Data collection complete. {len(all_data)} successful tests.")

    return pd.DataFrame(all_data)


def plot_cache_misses(df, iterations, output_file):
    """
    Create a plot showing cache misses vs array size.

    Args:
        df: DataFrame with cache statistics
        iterations: Filter data for this iteration count
        output_file: Output filename for the plot
    """
    # Filter data for specific iteration count
    data = df[df['iterations'] == iterations].copy()
    data = data.sort_values('array_size')

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(f'Cache Performance Analysis (Iterations: {iterations})',
                 fontsize=16, fontweight='bold')

    # Plot 1: Cache misses vs array size (log scale)
    ax1.plot(data['array_size'], data['l1_misses'],
             'o-', linewidth=2, markersize=8, label='L1 Cache Misses (64KB)', color='#e74c3c')
    ax1.plot(data['array_size'], data['l2_misses'],
             '^-', linewidth=2, markersize=8, label='L2 Cache Misses (512KB)', color='#f39c12')
    ax1.plot(data['array_size'], data['l3_misses'],
             's-', linewidth=2, markersize=8, label='L3 Cache Misses (256MB)', color='#3498db')

    ax1.set_xlabel('Array Size (elements)', fontsize=12)
    ax1.set_ylabel('Cache Misses', fontsize=12)
    ax1.set_title('Cache Misses vs Array Size', fontsize=14)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.legend(fontsize=10)

    # Add cache size reference lines (in elements, 1 double = 8 bytes)
    # L1: 64 KB = 8,192 elements
    # L2: 512 KB = 65,536 elements
    # L3: 256 MB = 33,554,432 elements
    ax1.axvline(x=8192, color='#e74c3c', linestyle='--', alpha=0.3, linewidth=1.5)
    ax1.axvline(x=65536, color='#f39c12', linestyle='--', alpha=0.3, linewidth=1.5)
    ax1.axvline(x=33554432, color='#3498db', linestyle='--', alpha=0.3, linewidth=1.5)

    # Plot 2: Miss rate (percentage)
    # Calculate miss rates using correct denominators for each cache level
    # L1 miss rate = l1_misses / l1_loads
    # L2 miss rate = l2_misses / l2_accesses
    # L3 miss rate = l3_misses / l3_accesses
    data['l1_miss_rate'] = ((data['l1_misses'] / data['l1_loads']) * 100).fillna(0)
    data['l2_miss_rate'] = ((data['l2_misses'] / data['l2_accesses']) * 100).fillna(0)
    data['l3_miss_rate'] = ((data['l3_misses'] / data['l3_accesses']) * 100).fillna(0)

    # Replace inf values with 0 (in case of division by zero)
    data['l1_miss_rate'] = data['l1_miss_rate'].replace([float('inf'), -float('inf')], 0)
    data['l2_miss_rate'] = data['l2_miss_rate'].replace([float('inf'), -float('inf')], 0)
    data['l3_miss_rate'] = data['l3_miss_rate'].replace([float('inf'), -float('inf')], 0)

    ax2.plot(data['array_size'], data['l1_miss_rate'],
             'o-', linewidth=2, markersize=8, label='L1 Miss Rate', color='#e74c3c')
    ax2.plot(data['array_size'], data['l2_miss_rate'],
             '^-', linewidth=2, markersize=8, label='L2 Miss Rate', color='#f39c12')
    ax2.plot(data['array_size'], data['l3_miss_rate'],
             's-', linewidth=2, markersize=8, label='L3 Miss Rate', color='#3498db')

    ax2.set_xlabel('Array Size (elements)', fontsize=12)
    ax2.set_ylabel('Miss Rate (%)', fontsize=12)
    ax2.set_title('Cache Miss Rate vs Array Size', fontsize=14)
    ax2.set_xscale('log')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved: {output_file}")
    plt.close()


def generate_report(df, output_file="cache_report.txt"):
    """
    Generate a text report with statistics.

    Args:
        df: DataFrame with cache statistics
        output_file: Output filename for the report
    """
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CACHE PERFORMANCE ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")

        for iterations in ITERATION_CONFIGS:
            data = df[df['iterations'] == iterations].copy()
            data = data.sort_values('array_size')

            f.write(f"\n{'='*80}\n")
            f.write(f"ITERATIONS: {iterations}\n")
            f.write(f"{'='*80}\n\n")

            f.write(f"{'Array Size':<15} {'Size (KB)':<12} {'L1 Misses':<15} {'L2 Misses':<15} {'L3 Misses':<15}\n")
            f.write("-" * 80 + "\n")

            for _, row in data.iterrows():
                f.write(f"{row['array_size']:<15,} {row['array_size_kb']:<12.1f} "
                       f"{row['l1_misses']:<15,} {row['l2_misses']:<15,} {row['l3_misses']:<15,}\n")

            f.write("\n")

    print(f"Report saved: {output_file}")


def main():
    """Main execution function."""
    print("=" * 60)
    print("Cache Miss Analysis Tool")
    print("=" * 60)
    print()

    # Check if executable exists
    if not os.path.exists(EXECUTABLE):
        print(f"Error: Executable '{EXECUTABLE}' not found!")
        print("Please build danial_cache.cpp first.")
        return 1

    # Check if perf is available
    try:
        subprocess.run(["perf", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: 'perf' command not found!")
        print("Please install linux-tools package.")
        return 1

    # Note about permissions
    print("Note: AMD L2 metrics may require elevated permissions.")
    print("If L2 counters are zero, try adjusting perf_event_paranoid:")
    print("  sudo sysctl kernel.perf_event_paranoid=1")
    print("L1 and L3 counters use generic events and should work without special permissions.")
    print()

    # Collect data
    df = collect_data()

    if df.empty:
        print("Error: No data collected!")
        return 1

    # Save raw data
    csv_file = "cache_data.csv"
    df.to_csv(csv_file, index=False)
    print(f"\nRaw data saved: {csv_file}")

    # Generate plots
    print("\nGenerating plots...")
    for iterations in ITERATION_CONFIGS:
        output_file = f"cache_misses_{iterations}iter.png"
        plot_cache_misses(df, iterations, output_file)

    # Generate report
    generate_report(df)

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
