import sys
import json
import csv
import time
import tracemalloc
from pathlib import Path
from datetime import datetime
import yaml
import fkcompute
import multiprocessing as mp
import click


def measure_fk_computation(knot, knot_config, degree):
    tracemalloc.start()
    start_time = time.time()

    try:
        # Get other parameters from the YAML
        braid = knot_config.get("braid")
        inversion = {"inversion_data": knot_config.get("inversion"), "braid": braid}
        threads = mp.cpu_count()
        save_data = True

        # Call fkcompute.fk() with individual parameters using the crossing number as degree
        fkcompute.fk(
            braid,
            degree=degree,
            name=knot,
            threads=threads,
            save_data=save_data,
            inversion=inversion,
        )
        success = True
        error = None
    except Exception as e:
        success = False
        error = str(e)

    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "success": success,
        "error": error,
        "time_seconds": end_time - start_time,
        "peak_memory_mb": peak / (1024 * 1024),
        "current_memory_mb": current / (1024 * 1024),
    }


def get_crossing_number(knot_name):
    """Extract crossing number from knot name (e.g., '10_100' -> 10, '11n_119' -> 11, '12a_1' -> 12)."""
    try:
        # Handle formats like '10_100', '11n_119', '12a_1', etc.
        prefix = knot_name.split('_')[0]
        # Extract leading digits from the prefix
        import re
        match = re.match(r'^(\d+)', prefix)
        if match:
            return int(match.group(1))
        return None
    except (ValueError, IndexError):
        return None


def run_benchmark(degrees, max_crossings=None, verbose=False, knots_file="all_knots.yaml"):
    results = []
    overall_start = time.time()

    with open(knots_file, "r") as f:
        all_knots = yaml.safe_load(f)

    # Filter knots by crossing number if specified
    if max_crossings is not None:
        filtered_knots = {
            k: v for k, v in all_knots.items()
            if get_crossing_number(k) is not None and get_crossing_number(k) <= max_crossings
        }
    else:
        filtered_knots = all_knots

    total_knots = len(filtered_knots)
    total_computations = total_knots * len(degrees)

    print(f"\n{'='*60}")
    print(f"FK Computation Run")
    print(f"{'='*60}")
    print(f"Knots: {total_knots} (max crossings: {max_crossings if max_crossings else 'all'})")
    print(f"Degrees: {degrees}")
    print(f"Total computations: {total_computations}")
    print(f"{'='*60}")

    for degree in degrees:
        print(f"\n{'='*60}")
        print(f"Computing degree {degree}")
        print(f"{'='*60}")

        for idx, knot in enumerate(filtered_knots, 1):
            knot_name = knot

            metrics = measure_fk_computation(knot, filtered_knots[knot], degree)

            result = {
                "degree": degree,
                "knot": knot_name,
                **metrics,
            }
            results.append(result)

            if verbose:
                if metrics["success"]:
                    print(
                        f"[{idx}/{total_knots}] {knot_name}: OK ({metrics['time_seconds']:.2f}s, {metrics['peak_memory_mb']:.1f}MB)"
                    )
                else:
                    print(f"[{idx}/{total_knots}] {knot_name}: FAILED: {metrics['error']}")
            else:
                status = "✓" if metrics["success"] else "✗"
                print(
                    f"  [{idx}/{total_knots}] {status} {knot_name:20s} deg={degree:2d} {metrics['time_seconds']:8.2f}s {metrics['peak_memory_mb']:10.1f}MB"
                )

    overall_end = time.time()
    total_time = overall_end - overall_start

    return results, total_time


def save_results_csv(results, filepath):
    if not results:
        return
    with open(filepath, "a", newline="") as f:
        fieldnames = [
            "degree",
            "knot",
            "success",
            "time_seconds",
            "peak_memory_mb",
            "current_memory_mb",
            "error",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writerows(results)

    print(f"\nResults saved to {filepath}")


def save_results_json(results, filepath):
    with open(filepath, "w") as f:
        json.dump(
            {"timestamp": datetime.now().isoformat(), "results": results}, f, indent=2
        )

    print(f"\nResults saved to {filepath}")


def format_time(seconds):
    """Format time in a human-readable way."""
    if seconds < 60:
        return f"{seconds:.2f}s"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.1f}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{hours}h {minutes}m {secs:.1f}s"


def print_summary(results, total_time):
    if not results:
        print("\nNo results to summarize")
        return

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")

    total_successful = sum(1 for r in results if r["success"])
    total_failed = len(results) - total_successful

    print(f"\nOverall Statistics:")
    print(f"  Total computations: {len(results)}")
    print(f"  Successful:         {total_successful}")
    print(f"  Failed:             {total_failed}")
    print(f"  Total time:         {format_time(total_time)}")
    print(f"  Average per knot:   {format_time(total_time / len(results))}")

    by_degree = {}
    for r in results:
        deg = r["degree"]
        if deg not in by_degree:
            by_degree[deg] = []
        by_degree[deg].append(r)

    for degree in sorted(by_degree.keys()):
        deg_results = by_degree[degree]
        successful = [r for r in deg_results if r["success"]]
        failed = [r for r in deg_results if not r["success"]]

        print(f"\n{'='*60}")
        print(f"Degree {degree}")
        print(f"{'='*60}")
        print(f"  Computations: {len(deg_results)} ({len(successful)} successful, {len(failed)} failed)")

        if successful:
            total_time_deg = sum(r["time_seconds"] for r in successful)
            avg_time = total_time_deg / len(successful)
            max_time = max(r["time_seconds"] for r in successful)
            min_time = min(r["time_seconds"] for r in successful)
            avg_memory = sum(r["peak_memory_mb"] for r in successful) / len(successful)
            max_memory = max(r["peak_memory_mb"] for r in successful)

            print(f"  Total time:   {format_time(total_time_deg)}")
            print(f"  Avg time:     {format_time(avg_time)}")
            print(f"  Min time:     {format_time(min_time)}")
            print(f"  Max time:     {format_time(max_time)}")
            print(f"  Avg memory:   {avg_memory:8.1f}MB")
            print(f"  Max memory:   {max_memory:8.1f}MB")

            if failed:
                print(f"\n  Failed knots:")
                for r in failed:
                    print(f"    - {r['knot']}: {r['error']}")


@click.command()
@click.option(
    '--max-crossings',
    type=int,
    required=True,
    help='Maximum number of crossings for knots to compute (e.g., 10 for knots up to 10 crossings)'
)
@click.option(
    '--degree',
    type=int,
    required=True,
    multiple=True,
    help='Degree(s) to compute FK to. Can be specified multiple times (e.g., --degree 10 --degree 12)'
)
@click.option(
    '--output',
    type=click.Path(),
    help='Output file path for results (CSV format). If not specified, results will only be printed.'
)
@click.option(
    '--json-output',
    type=click.Path(),
    help='Output file path for JSON results. If not specified, only CSV/print output will be used.'
)
@click.option(
    '--verbose',
    is_flag=True,
    help='Enable verbose output with detailed information for each computation'
)
@click.option(
    '--knots-file',
    type=click.Path(exists=True),
    default='all_knots.yaml',
    help='Path to YAML file containing knot configurations (default: all_knots.yaml)'
)
def main(max_crossings, degree, output, json_output, verbose, knots_file):
    """
    Compute FK homology for knots up to a certain number of crossings.

    This tool computes FK homology for knots and tracks detailed statistics
    about computation time and memory usage for each knot and overall.

    Examples:

      # Compute FK for knots up to 10 crossings to degree 10
      python compute_fks.py --max-crossings 10 --degree 10

      # Compute for multiple degrees and save results
      python compute_fks.py --max-crossings 8 --degree 8 --degree 10 --output results.csv

      # Verbose output with JSON export
      python compute_fks.py --max-crossings 12 --degree 12 --verbose --json-output results.json
    """
    degrees = list(degree)

    if not degrees:
        click.echo("Error: At least one --degree must be specified", err=True)
        sys.exit(1)

    click.echo(f"Starting FK computations")
    click.echo(f"Max crossings: {max_crossings}")
    click.echo(f"Degrees: {degrees}")

    results, total_time = run_benchmark(
        degrees=degrees,
        max_crossings=max_crossings,
        verbose=verbose,
        knots_file=knots_file
    )

    print_summary(results, total_time)

    if output:
        save_results_csv(results, output)

    if json_output:
        save_results_json(results, json_output)


if __name__ == "__main__":
    main()
