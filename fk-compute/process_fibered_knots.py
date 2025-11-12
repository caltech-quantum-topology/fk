#!/usr/bin/env python3

import json
import os
import sys
import traceback
import fkcompute
from pathlib import Path

def load_fibered_table(filepath="fibered_table.json"):
    """Load knot data from fibered_table.json"""
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: Could not find {filepath}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in {filepath}: {e}")
        sys.exit(1)

def parse_braid_string(braid_str):
    """Convert braid string like '{1,1,-2,1,3,2,2,-4,-3,2,-3,-4}' to list"""
    try:
        braid_str = braid_str.strip('{}')
        return [int(x.strip()) for x in braid_str.split(',')]
    except ValueError as e:
        raise ValueError(f"Failed to parse braid string '{braid_str}': {e}")

def compute_fk_for_knot(knot_data, results_dir="results"):
    """Compute FK invariant for a single knot and save result"""
    knot_name = knot_data["name"]
    braid_str = knot_data["braid"]
    max_x_degree = int(knot_data["max_x_degree"])+3

    print(f"Processing knot {knot_name}...")

    try:
        # Parse braid string to list
        braid_list = parse_braid_string(braid_str)

        # Import and use the fkcompute module
        from fkcompute import fk

        # Compute FK invariant
        result = fk(braid_list, max_x_degree, verbose=False)

        # Save result to file
        result_file = os.path.join(results_dir, f"{knot_name}.json")
        with open(result_file, 'w') as f:
            json.dump(result, f, indent=2)

        print(f"✓ Successfully computed and saved {knot_name}")
        return True

    except Exception as e:
        error_msg = f"✗ Failed to compute {knot_name}: {str(e)}"
        print(error_msg)

        # Save error to a separate file for debugging
        error_file = os.path.join(results_dir, f"{knot_name}_error.txt")
        with open(error_file, 'w') as f:
            f.write(f"Error computing FK for {knot_name}:\n")
            f.write(f"Braid: {braid_str}\n")
            f.write(f"Max degree: {max_x_degree}\n")
            f.write(f"Error: {str(e)}\n\n")
            f.write("Full traceback:\n")
            f.write(traceback.format_exc())

        return False

def main():
    """Main function to process all knots in fibered_table.json"""
    # Ensure results directory exists
    results_dir = "results"
    Path(results_dir).mkdir(exist_ok=True)

    # Load knot data
    print("Loading fibered knot data...")
    knots = load_fibered_table()
    print(f"Found {len(knots)} knots to process")

    # Process each knot
    successful = 0
    failed = 0

    for i, knot in enumerate(knots, 1):
        print(f"\n[{i}/{len(knots)}] ", end="")

        if compute_fk_for_knot(knot, results_dir):
            successful += 1
        else:
            failed += 1

    # Summary
    print(f"\n" + "="*50)
    print(f"Processing complete!")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Total: {len(knots)}")
    print(f"Results saved in: {results_dir}/")

    if failed > 0:
        print(f"Error details saved in *_error.txt files")

if __name__ == "__main__":
    main()
