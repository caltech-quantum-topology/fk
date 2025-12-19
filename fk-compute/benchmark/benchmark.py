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


def run_benchmark(degrees, file_pattern=None, verbose=False):
    results = []

    for degree in degrees:

        with open("all_knots.yaml", "r") as f:
            all_knots = yaml.safe_load(f)

        print(f"\n{'='*60}")
        print(f"Degree {degree}: {len(all_knots)} configurations")
        print(f"{'='*60}")

        for knot in all_knots:
            knot_name = knot

            metrics = measure_fk_computation(knot, all_knots[knot], degree)

            result = {
                "degree": degree,
                "knot": knot_name,
                **metrics,
            }
            results.append(result)

            if verbose:
                if metrics["success"]:
                    print(
                        f"OK ({metrics['time_seconds']:.2f}s, {metrics['peak_memory_mb']:.1f}MB)"
                    )
                else:
                    print(f"FAILED: {metrics['error']}")
            else:
                status = "✓" if metrics["success"] else "✗"
                print(
                    f"  {status} {knot_name:20s} deg={degree:2d} {metrics['time_seconds']:8.2f}s {metrics['peak_memory_mb']:10.1f}MB"
                )

    return results


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


def print_summary(results):
    if not results:
        print("\nNo results to summarize")
        return

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")

    by_degree = {}
    for r in results:
        deg = r["degree"]
        if deg not in by_degree:
            by_degree[deg] = []
        by_degree[deg].append(r)

    for degree in sorted(by_degree.keys()):
        deg_results = by_degree[degree]
        successful = [r for r in deg_results if r["success"]]

        if successful:
            total_time = sum(r["time_seconds"] for r in successful)
            avg_time = total_time / len(successful)
            max_time = max(r["time_seconds"] for r in successful)
            avg_memory = sum(r["peak_memory_mb"] for r in successful) / len(successful)
            max_memory = max(r["peak_memory_mb"] for r in successful)

            print(f"\nDegree {degree}: {len(successful)}/{len(deg_results)} successful")
            print(f"  Total time:   {total_time:8.2f}s")
            print(f"  Avg time:     {avg_time:8.2f}s")
            print(f"  Max time:     {max_time:8.2f}s")
            print(f"  Avg memory:   {avg_memory:8.1f}MB")
            print(f"  Max memory:   {max_memory:8.1f}MB")


def main():
    degrees = [10]
    output = None #"benchmark_results.csv"
    print(f"Benchmarking FK computations for degrees: {degrees}")

    results = run_benchmark(
        degrees=degrees,
    )

    print_summary(results)

    if output:
        save_results_csv(results, output)
    if results:
        fieldnames = ["degree", "knot", "success", "time_seconds", "peak_memory_mb"]
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
