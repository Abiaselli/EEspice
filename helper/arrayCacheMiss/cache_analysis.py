#!/home/boris/github/EEspice/venv/bin/python3
"""
Cache Miss Analysis Script

This script runs danial_cache.o with varying array sizes and collects
L1, L2, and L3 cache miss statistics using perf. It then generates plots
showing the relationship between array size and cache misses.

Supports both Intel and AMD CPUs with appropriate performance counters:
- Intel: Uses l2_rqsts events for L2 cache statistics
- AMD: Uses all_l2_cache_* metrics with -M flag
- Generic: L1-dcache and LLC (L3) events work on both platforms

L3 miss rate is calculated as cache-misses / cache-references.

Usage:
    python cache_analysis.py                    # Auto-detect CPU
    python cache_analysis.py --cpu-vendor intel # Force Intel configuration
    python cache_analysis.py --verbose          # Show detailed perf events
    python cache_analysis.py --use-generic      # Use only generic events
"""

import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import argparse

# Configuration
EXECUTABLE = "./danial_cache.o"
ARRAY_SIZES = [
    1000, 2000, 4000, 8000,           # < 64KB (fits in L0-D)
    16000, 30000, 50000,              # ~128KB-400KB (exceeds L0-D+L1, fits in L2)
    100000, 200000, 400000,           # ~800KB-3.2MB (around L2 boundary)
    600000, 1000000, 2000000,         # ~4.8MB-16MB (exceeds L2, fits in L3)
    5000000, 10000000                 # > 40MB (exceeds L3)
]
ITERATION_CONFIGS = [1, 100000]  # Test with 1 iteration and 100000 iterations

# Global configuration
CPU_VENDOR = "auto"  # Will be set by detect_cpu_vendor() or command-line
VERBOSE = False      # Will be set by command-line


def detect_cpu_vendor():
    """
    Detect CPU vendor from /proc/cpuinfo.

    Returns:
        str: 'intel', 'amd', or 'unknown'
    """
    try:
        with open('/proc/cpuinfo', 'r') as f:
            cpuinfo = f.read().lower()

        if 'genuineintel' in cpuinfo:
            return 'intel'
        elif 'authenticamd' in cpuinfo:
            return 'amd'
        else:
            # Check vendor_id line as fallback
            for line in cpuinfo.split('\n'):
                if line.startswith('vendor_id'):
                    if 'intel' in line:
                        return 'intel'
                    elif 'amd' in line:
                        return 'amd'
            return 'unknown'
    except Exception as e:
        if VERBOSE:
            print(f"Warning: Could not detect CPU vendor: {e}")
        return 'unknown'


def get_perf_events(cpu_vendor):
    """
    Get appropriate perf events based on CPU vendor.

    Args:
        cpu_vendor: 'intel', 'amd', or 'unknown'

    Returns:
        tuple: (events_string, metrics_string, use_metrics_flag)
    """
    # L1 and L3 events are generic and work on both Intel and AMD
    base_events = "L1-dcache-loads,L1-dcache-load-misses,cache-references,cache-misses"

    if cpu_vendor == 'intel':
        # For Intel, we'll use generic events since specific L2 events may not be available
        # Note: On modern Intel systems, cache-references/misses typically refer to LLC (L3)
        # L2 statistics may not be directly accessible via perf on all Intel CPUs
        events = base_events
        metrics = ""  # Intel doesn't use -M flag
        use_metrics = False

        if VERBOSE:
            print(f"Using Intel performance events:")
            print(f"  Events: {base_events}")
            print(f"  Note: L2 statistics estimated from L1 and L3 miss patterns")

    elif cpu_vendor == 'amd':
        # AMD-specific configuration
        events = base_events
        metrics = "all_l2_cache_accesses,all_l2_cache_misses"
        use_metrics = True

        if VERBOSE:
            print(f"Using AMD performance events:")
            print(f"  Base events: {base_events}")
            print(f"  L2 metrics: {metrics}")

    else:
        # Unknown vendor - use generic events only
        events = base_events
        metrics = ""
        use_metrics = False

        if VERBOSE:
            print(f"Using generic performance events only (unknown CPU vendor)")
            print(f"  Events: {base_events}")

    return events, metrics, use_metrics


def run_perf_test(array_size, iterations, cpu_vendor):
    """
    Run danial_cache.o with perf stat to collect cache statistics.

    Args:
        array_size: Size of array to test
        iterations: Number of iterations to run
        cpu_vendor: CPU vendor ('intel', 'amd', or 'unknown')

    Returns:
        Dictionary with cache statistics or None if failed
    """
    # Get appropriate perf events for this CPU
    events, metrics, use_metrics = get_perf_events(cpu_vendor)

    # Build perf command
    cmd = [
        "perf", "stat",
        "--no-big-num",            # Disable locale/group separators
        "-x", ";",                 # CSV-like output with semicolon delimiter
        "-e", events
    ]

    # Only add -M flag for AMD CPUs with metrics
    if use_metrics and metrics:
        cmd.extend(["-M", metrics])

    cmd.extend([
        EXECUTABLE,
        str(array_size),
        str(iterations)
    ])

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
        stats = parse_perf_output(result.stderr, cpu_vendor)
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


def parse_perf_output(output: str, cpu_vendor: str):
    """
    Parse perf stat CSV output to extract cache statistics.

    With --no-big-num and -x ";", perf stat outputs CSV format:
    - Event lines: value;unit;event;runtime;pcnt;[var];[metric_value];[metric_unit]
    - Metric lines: Additional metrics with empty earlier fields, only metric_value;metric_name filled

    Args:
        output: stderr output from perf stat (CSV format)
        cpu_vendor: CPU vendor ('intel', 'amd', or 'unknown')

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

    # For Intel CPUs without specific L2 events, estimate L2 behavior
    # L2 accesses ≈ L1 misses (what misses L1 goes to L2)
    # L2 misses ≈ L3 accesses (what misses L2 goes to L3)
    if cpu_vendor == 'intel' and stats['l2_accesses'] == 0 and stats['l2_misses'] == 0:
        stats['l2_accesses'] = stats['l1_misses']
        stats['l2_misses'] = stats['l3_accesses']

    return stats


def collect_data(cpu_vendor):
    """
    Collect cache miss data for all array sizes and iteration counts.

    Args:
        cpu_vendor: CPU vendor ('intel', 'amd', or 'unknown')

    Returns:
        pandas DataFrame with all collected data
    """
    all_data = []
    total_tests = len(ARRAY_SIZES) * len(ITERATION_CONFIGS)
    test_count = 0

    print(f"Running {total_tests} tests on {cpu_vendor.upper()} CPU...")
    print(f"Array sizes: {len(ARRAY_SIZES)} configurations")
    print(f"Iteration counts: {ITERATION_CONFIGS}")
    print("-" * 60)

    for iterations in ITERATION_CONFIGS:
        print(f"\nTesting with {iterations} iteration(s):")
        for array_size in ARRAY_SIZES:
            test_count += 1
            print(f"  [{test_count}/{total_tests}] Array size: {array_size:>8} elements "
                  f"({(array_size*8/1024):.1f} KB)...", end=" ", flush=True)

            stats = run_perf_test(array_size, iterations, cpu_vendor)

            if stats:
                all_data.append(stats)
                print(f"L1: {stats['l1_misses']:>10,}, L2: {stats['l2_misses']:>10,}, L3: {stats['l3_misses']:>10,}")
            else:
                print("FAILED")

    print("\n" + "=" * 60)
    print(f"Data collection complete. {len(all_data)} successful tests.")

    return pd.DataFrame(all_data)


def plot_cache_misses(df, iterations, output_file, cpu_vendor):
    """
    Create a plot showing cache misses vs array size.

    Args:
        df: DataFrame with cache statistics
        iterations: Filter data for this iteration count
        output_file: Output filename for the plot
        cpu_vendor: CPU vendor for cache size labels
    """
    # Filter data for specific iteration count
    data = df[df['iterations'] == iterations].copy()
    data = data.sort_values('array_size')

    # Set cache size labels based on CPU vendor
    if cpu_vendor == 'intel':
        # Intel cache sizes for modern CPUs
        # L0-D (48KB) + L1 (192KB) = 240KB total per-core data cache
        l1_label = 'L1 Data Cache Misses (240KB per core)'
        l2_label = 'L2 Cache Misses [estimated] (3MB per core)'
        l3_label = 'L3 Cache Misses (36MB shared)'
        # Intel cache sizes in elements (doubles = 8 bytes)
        l1_size = 30720     # 240KB / 8 bytes (L0-D + L1 data cache)
        l2_size = 393216    # 3MB / 8 bytes
        l3_size = 4718592   # 36MB / 8 bytes
    else:
        # AMD cache sizes (original values)
        l1_label = 'L1 Cache Misses (64KB)'
        l2_label = 'L2 Cache Misses (512KB)'
        l3_label = 'L3 Cache Misses (256MB)'
        l1_size = 8192      # 64KB / 8 bytes
        l2_size = 65536     # 512KB / 8 bytes
        l3_size = 33554432  # 256MB / 8 bytes

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(f'Cache Performance Analysis ({cpu_vendor.upper()} CPU, Iterations: {iterations})',
                 fontsize=16, fontweight='bold')

    # Plot 1: Cache misses vs array size (log scale)
    ax1.plot(data['array_size'], data['l1_misses'],
             'o-', linewidth=2, markersize=8, label=l1_label, color='#e74c3c')
    ax1.plot(data['array_size'], data['l2_misses'],
             '^-', linewidth=2, markersize=8, label=l2_label, color='#f39c12')
    ax1.plot(data['array_size'], data['l3_misses'],
             's-', linewidth=2, markersize=8, label=l3_label, color='#3498db')

    ax1.set_xlabel('Array Size (elements)', fontsize=12)
    ax1.set_ylabel('Cache Misses', fontsize=12)
    ax1.set_title('Cache Misses vs Array Size', fontsize=14)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.legend(fontsize=10)

    # Add cache size reference lines (in elements, 1 double = 8 bytes)
    ax1.axvline(x=l1_size, color='#e74c3c', linestyle='--', alpha=0.3, linewidth=1.5)
    ax1.axvline(x=l2_size, color='#f39c12', linestyle='--', alpha=0.3, linewidth=1.5)
    ax1.axvline(x=l3_size, color='#3498db', linestyle='--', alpha=0.3, linewidth=1.5)

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
    global CPU_VENDOR, VERBOSE

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Cache miss analysis tool for Intel and AMD CPUs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                    # Auto-detect CPU and run analysis
  %(prog)s --cpu-vendor intel # Force Intel CPU configuration
  %(prog)s --verbose         # Show detailed perf event information
  %(prog)s --use-generic     # Use only generic perf events
        """
    )
    parser.add_argument('--cpu-vendor', choices=['auto', 'intel', 'amd'],
                        default='auto',
                        help='CPU vendor (default: auto-detect)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output for debugging')
    parser.add_argument('--use-generic', action='store_true',
                        help='Use only generic perf events (no vendor-specific events)')

    args = parser.parse_args()

    VERBOSE = args.verbose

    print("=" * 60)
    print("Cache Miss Analysis Tool")
    print("=" * 60)
    print()

    # Check if executable exists
    if not os.path.exists(EXECUTABLE):
        print(f"Error: Executable '{EXECUTABLE}' not found!")
        print("Please build danial_cache.cpp first.")
        print("Build command: g++ -O2 danial_cache.cpp -o danial_cache.o")
        return 1

    # Check if perf is available
    try:
        subprocess.run(["perf", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: 'perf' command not found!")
        print("Please install linux-tools package.")
        print("  Ubuntu/Debian: sudo apt-get install linux-tools-common linux-tools-generic")
        print("  Fedora/RHEL: sudo dnf install perf")
        return 1

    # Determine CPU vendor
    if args.use_generic:
        cpu_vendor = 'unknown'  # Forces generic events only
        print("Using generic perf events only (--use-generic flag)")
    elif args.cpu_vendor == 'auto':
        cpu_vendor = detect_cpu_vendor()
        print(f"Detected CPU vendor: {cpu_vendor.upper()}")
        if cpu_vendor == 'unknown':
            print("Warning: Could not detect CPU vendor, using generic events only")
    else:
        cpu_vendor = args.cpu_vendor
        print(f"Using specified CPU vendor: {cpu_vendor.upper()}")

    # Note about permissions
    if cpu_vendor == 'amd':
        print("\nNote: AMD L2 metrics may require elevated permissions.")
        print("If L2 counters are zero, try adjusting perf_event_paranoid:")
        print("  sudo sysctl kernel.perf_event_paranoid=1")
    elif cpu_vendor == 'intel':
        print("\nNote: For Intel CPUs, L2 statistics are estimated from L1/L3 miss patterns")
        print("as direct L2 events may not be available on all Intel systems.")
        print("L2 accesses ≈ L1 misses, L2 misses ≈ L3 accesses")
    print("L1 and L3 counters use generic events and should work without special permissions.")
    print()

    # Collect data
    df = collect_data(cpu_vendor)

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
        output_file = f"cache_misses_{cpu_vendor}_{iterations}iter.png"
        plot_cache_misses(df, iterations, output_file, cpu_vendor)

    # Generate report
    generate_report(df)

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
