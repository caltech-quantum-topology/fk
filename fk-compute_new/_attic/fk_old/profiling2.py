from memory_profiler import memory_usage
from main import *
import json
import random
import sys
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import threading
import queue
import time

def load_knot_data():
    """Load knot data from inversion_data.json and filter by braid length 3-12"""
    with open('Data/Input/inversion_data.json', 'r') as f:
        data = json.load(f)

    # Group knots by braid length
    knots_by_length = {}
    for knot in data:
        if 'braid' in knot and knot['braid']:
            braid_length = len(knot['braid'])
            if 3 <= braid_length <= 12:
                if braid_length not in knots_by_length:
                    knots_by_length[braid_length] = []
                knots_by_length[braid_length].append(knot['braid'])

    return knots_by_length

def increase_length(braid):
    return braid + [2, -2]

def get_sign_assignment_wrapper(braid, degree, assignment, verbose):
    assignment.update(get_sign_assignment(braid, degree, verbose))

def compute_wrapper(braid_states, degree, name):
    braid_ilp.save(braid_states=braid_states, degree=degree, save_to=f'Data/Input/{name}.csv')
    os.system(f'./_ "Data/Input/{name}" "Data/Output/{name}"')

def process_single_knot(args):
    """Process a single knot across all degrees"""
    braid, knot_index, braid_length, degrees = args

    results = []

    # Perform sign assignment search for this braid
    start_time = time.time()
    assignment = {}
    try:
        mem_usage = memory_usage((lambda: get_sign_assignment_wrapper(braid, degree=10, assignment=assignment, verbose=False)), interval=0.1)
        end_time = time.time()
        sign_time = end_time - start_time
        sign_max_mem = max(mem_usage)

        if 'braid' not in assignment or 'inversion_data' not in assignment:
            return []

        # Set up braid states with sign assignment
        braid = assignment['braid']
        bs = BraidStates(braid)
        load_sign_data(bs, assignment['inversion_data'], compute_r_matrices=True)

        # Test across all degrees for this braid
        for degree in degrees:
            start_time = time.time()
            name = f"braid_len_{braid_length}_deg_{degree}_knot_{knot_index}"
            mem_usage = memory_usage((lambda: compute_wrapper(braid_states=bs, degree=degree, name=name)), interval=0.1)
            end_time = time.time()
            comp_time = end_time - start_time
            comp_max_mem = max(mem_usage)

            results.append((braid_length, degree, sign_time, sign_max_mem, comp_time, comp_max_mem, knot_index))
    except Exception as e:
        print(f"Error processing knot {knot_index} of length {braid_length}: {e}")
        return []

    return results

if __name__ == "__main__":
    import time

    # Parse command line arguments for multiple file option
    file_suffix = ""
    num_workers = mp.cpu_count()  # Use all available cores
    knots_per_length = 3  # Process multiple knots per length

    if len(sys.argv) > 1:
        file_suffix = f"_{sys.argv[1]}"
    if len(sys.argv) > 2:
        num_workers = int(sys.argv[2])
    if len(sys.argv) > 3:
        knots_per_length = int(sys.argv[3])

    print(f"Using {num_workers} workers, processing {knots_per_length} knots per length")

    # Load knot data
    print("Loading knot data...")
    knots_by_length = load_knot_data()

    degrees = list(range(5, 33))
    #degrees = list(range(5, 6))

    output_filename = f"profiling_results2{file_suffix}.txt"

    # Prepare all tasks
    tasks = []
    for braid_length in range(3, 13):  # Lengths 3 through 12
        if braid_length not in knots_by_length:
            print(f"No knots found for length {braid_length}, skipping...")
            continue

        available_braids = knots_by_length[braid_length]
        print(f"Found {len(available_braids)} knots of length {braid_length}")

        # Sample multiple braids for this length
        num_to_sample = min(knots_per_length, len(available_braids))
        selected_braids = random.sample(available_braids, num_to_sample)

        for i, braid in enumerate(selected_braids):
            knot_index = available_braids.index(braid)
            tasks.append((braid, knot_index, braid_length, degrees))

    print(f"Processing {len(tasks)} knots in parallel...")

    # Thread-safe file writing
    write_lock = threading.Lock()
    results_queue = queue.Queue()

    def write_results():
        with open(output_filename, "w") as f:
            f.write("Braid_Length,Degree,Sign_Assignment_Time,Sign_Assignment_Memory,Computation_Time,Computation_Memory,Knot_Index\n")

            while True:
                try:
                    result = results_queue.get(timeout=1)
                    if result is None:  # Sentinel to stop
                        break

                    braid_length, degree, sign_time, sign_max_mem, comp_time, comp_max_mem, knot_index = result
                    f.write(f"{braid_length},{degree},{sign_time:.2f},{sign_max_mem:.2f},{comp_time:.2f},{comp_max_mem:.2f},{knot_index}\n")
                    f.flush()
                except queue.Empty:
                    continue

    # Start the writer thread
    writer_thread = threading.Thread(target=write_results)
    writer_thread.start()

    # Process tasks in parallel
    completed_tasks = 0
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_task = {executor.submit(process_single_knot, task): task for task in tasks}

        for future in as_completed(future_to_task):
            task = future_to_task[future]
            braid, knot_index, braid_length, degrees = task

            try:
                results = future.result()
                for result in results:
                    results_queue.put(result)

                completed_tasks += 1
                print(f"Completed {completed_tasks}/{len(tasks)}: knot {knot_index} of length {braid_length}")

            except Exception as e:
                print(f"Error processing knot {knot_index} of length {braid_length}: {e}")

    # Signal writer thread to stop
    results_queue.put(None)
    writer_thread.join()

    print(f"Profiling completed. Results saved to {output_filename}")
