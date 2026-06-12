from memory_profiler import memory_usage
from main import *

def increase_length(braid):
    return braid + [2, -2]

def get_sign_assignment_wrapper(braid, degree, assignment, verbose):
    assignment.update(get_sign_assignment(braid, degree, verbose))

def compute_wrapper(braid_states, degree, name):
    braid_ilp.save(braid_states=braid_states, degree=degree, save_to=f'Data/Input/{name}.csv')
    os.system(f'./_ "Data/Input/{name}" "Data/Output/{name}"')

if __name__ == "__main__":
    braid = [1,1,1,-2,-1,-1,-1,-2]
    name = "8_20"

    degrees = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
    n_trials = 4

    import time

    with open("profiling_results2.txt", "w") as f:
        for _ in range(n_trials):
            start_time = time.time()
            assignment = {}
            mem_usage = memory_usage((lambda: get_sign_assignment_wrapper(braid, degree=10, assignment=assignment, verbose=False)), interval=0.1)
            end_time = time.time()
            max_mem = max(mem_usage)
            f.write(f"Time taken for sign assignment at length {len(braid)}: {end_time - start_time:.2f} seconds, Max memory: {max_mem} MiB\n")
            braid = assignment['braid']
            bs = BraidStates(braid)
            load_sign_data(bs, assignment['inversion_data'], compute_r_matrices=True)
            for degree in degrees:
                print(f"Profiling with braid length {len(braid)} and degree {degree}")
                start_time = time.time()
                mem_usage = memory_usage((lambda: compute_wrapper(braid_states=bs, degree=degree, name=name)), interval=0.1)
                end_time = time.time()
                max_mem = max(mem_usage)
                f.write(f"Braid length: {len(braid)}, Degree: {degree}, Time taken: {end_time - start_time:.2f} seconds, Max memory: {max_mem} MiB\n")
            braid = increase_length(braid) # Increase braid length for next iteration
