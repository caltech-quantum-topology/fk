from memory_profiler import memory_usage
from main import *
import json
import random
import sys
import time

def load_knot_data():
    """Load knot data from topology.fyi and filter by braid length 3-12"""
    data = load_from_topology_fyi()

    # Group knots by braid length
    knots_by_length = {}
    for knot in data:
        if 'braid' in knot and knot['braid']:
            braid_length = len(knot['braid'])
            if 3 <= braid_length <= 12:
                if braid_length not in knots_by_length:
                    knots_by_length[braid_length] = []
                knots_by_length[braid_length].append(knot)

    return knots_by_length

def load_from_topology_fyi ():
    import requests

    url = "https://topology.fyi/api/fk_fibered_inversion_data"

    headers = {
        "Accept": "application/json",
        "Authorization": "Basic dG9wb2xvZ3k6Znlp"
    }

    response = requests.get(url, headers=headers)
    data = response.json()

    return data

def increase_length(braid):
    return braid + [2, -2]

def compute_wrapper(braid_states, degree, name):
    braid_ilp.save(braid_states=braid_states, degree=degree, save_to=f'Data/Input/{name}.csv')
    os.system(f'./_ "Data/Input/{name}" "Data/Output/{name}"')

def process_single_knot_degree(args):
    """Process a single knot for a single degree"""
    knot_data, knotinfo_id, braid_length, degree = args

    try:
        # Use pre-computed sign assignment from topology data
        braid = knot_data['braid']
        inversion_data = knot_data['sign_assignment']

        # Set up braid states with pre-computed sign assignment
        bs = BraidStates(braid)
        load_sign_data(bs, inversion_data, compute_r_matrices=True)

        # Test for this specific degree
        start_time = time.time()
        name = f"braid_len_{braid_length}_deg_{degree}_knot_{knotinfo_id}"
        mem_usage = memory_usage((lambda: compute_wrapper(braid_states=bs, degree=degree, name=name)), interval=0.1)
        end_time = time.time()
        comp_time = end_time - start_time
        comp_max_mem = max(mem_usage)

        return (braid_length, degree, comp_time, comp_max_mem, knotinfo_id)
    except Exception as e:
        print(f"Error processing knot {knotinfo_id} of length {braid_length} at degree {degree}: {e}")
        return None

if __name__ == "__main__":
    import time
    # Parse command line arguments for multiple file option
    file_suffix = ""
    knots_per_length = 3  # Process multiple knots per length

    if len(sys.argv) > 1:
        file_suffix = f"_{sys.argv[1]}"
    if len(sys.argv) > 2:
        knots_per_length = int(sys.argv[2])

    print(f"Processing {knots_per_length} knots per length sequentially")

    # Load knot data
    print("Loading knot data...")
    knots_by_length = load_knot_data()

    degrees = list(range(5, 27, 3))

    output_filename = f"profiling_results_3{file_suffix}.txt"

    # Open output file for writing
    with open(output_filename, "w") as f:
        f.write("Braid_Length,Degree,Computation_Time_seconds,Computation_Memory_MB,KnotInfo_ID\n")

        # Process by degree first, then sample knots for each degree
        total_tasks = 0
        completed_tasks = 0

        for degree in degrees:
            print(f"Processing degree {degree}...")

            # Sample knots for this degree across all braid lengths
            degree_tasks = []
            for braid_length in range(3, 13):  # Lengths 3 through 12
                if braid_length not in knots_by_length:
                    continue

                available_knots = knots_by_length[braid_length]

                # Sample knots for this specific degree and braid length
                num_to_sample = min(knots_per_length, len(available_knots))
                selected_knots = random.sample(available_knots, num_to_sample)

                for knot_data in selected_knots:
                    knotinfo_id = knot_data.get('knotinfo_id', 'unknown')
                    degree_tasks.append((knot_data, knotinfo_id, braid_length, degree))
                    total_tasks += 1

            print(f"Processing {len(degree_tasks)} knots for degree {degree}")

            # Process all knots for this degree
            for task in degree_tasks:
                knot_data, knotinfo_id, braid_length, degree = task

                result = process_single_knot_degree(task)
                if result is not None:
                    braid_length, degree, comp_time, comp_max_mem, knotinfo_id = result
                    f.write(f"{braid_length},{degree},{comp_time:.2f},{comp_max_mem:.2f},{knotinfo_id}\n")
                    f.flush()

                completed_tasks += 1
                print(f"Completed {completed_tasks}: knot {knotinfo_id} of length {braid_length} at degree {degree}")

    print(f"Profiling completed. Results saved to {output_filename}")
