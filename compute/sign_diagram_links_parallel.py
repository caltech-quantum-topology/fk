from braidstates_links import BraidStates
from relations_links import full_reduce, check_sign_assignment
from braids import is_homogeneous_braid
import argparse
import os



BRAID_TYPE_HOMEGENOUS = 0
BRAID_TYPE_FIBERED = 1

def braid_type(braid):
    if is_homogeneous_braid(braid):
        return BRAID_TYPE_HOMEGENOUS
    else:
        return BRAID_TYPE_FIBERED

def flip(braid):
    n = BraidStates(braid).n_strands
    return list(map(lambda x: (2 * (x > 0) - 1) * (n - abs(x)), braid))


def rot(braid):
    """Rotate braid left by one position."""
    return braid[1:] + [braid[0]]

def generate_braid_variants(braid, include_flip=True, max_shifts=None):
    """
    Yield braid variants in order:
      shift=0: braid, braid flipped
      shift=1: rotated once, rotated once flipped
      ...
    Stops after len(braid) shifts, or earlier if max_shifts is set.
    """
    limit = len(braid) if max_shifts is None else min(len(braid), max_shifts)
    
    for shift in range(limit):
        yield braid, {"flipped": False, "shift": shift}
        if include_flip:
            yield flip(braid), {"flipped": True, "shift": shift}
        braid = rot(braid)



# 2) total_assignments_for_braid
def total_assignments_for_braid(braid_states):
    """
    Return total number of sign assignments (2^n_s_total) for this braid.
    """
    return 1 << braid_states.n_s_total


# 3) split_ranges
def split_ranges(n, chunk_size):
    """
    Yield half-open integer ranges [start, end) covering 0..n-1.
    """
    start = 0
    while start < n:
        end = min(start + chunk_size, n)
        yield start, end
        start = end

# 4) sign_assignment_from_index
def sign_assignment_from_index(index, braid_states):
    """
    Convert an integer index to a list of per-component sign arrays, matching:
        2*bit - 1  in the same order as in your current code.
    Returns a *new* list of lists (component-wise signs).
    """
    bits = list(map(int, list(bin(index)[2:].zfill(braid_states.n_s_total))))
    signs = list(map(lambda x: 2 * x - 1, bits))

    out = []
    cursor = 0
    for component in range(braid_states.n_components):
        k = braid_states.n_s[component]
        out.append(signs[cursor:cursor+k])
        cursor += k
    return out

# 5) check_assignment_for_braid
def check_assignment_for_braid(index, degree, braid, braid_states):
    """
    Try a single sign-assignment index for 'braid'.
    If valid, return (braid, strand_signs); else return None.
    """

    # fill strand signs from index
    per_component = sign_assignment_from_index(index, braid_states)
    for component in range(braid_states.n_components):
        braid_states.strand_signs[component] = per_component[component]

    # compute & validate
    braid_states.compute_matrices()
    if not braid_states.validate():
        return None

    braid_states.generate_position_assignments()
    all_relations = braid_states.get_state_relations()
    relations = full_reduce(all_relations)
    out = check_sign_assignment(degree, relations, braid_states)

    if out is not None:
        # Return the *exact* sign structure used (deep copy not strictly necessary for your flow)
        return braid, braid_states.strand_signs
    return None

# 6) worker_on_range
def worker_on_range(args):
    """
    Worker for a *range* of indices. Returns the first success in [start, end) or None.
    Args is a tuple: (degree, braid, start, end, progress)
      - 'progress' is a bool. If True, prints progress percentage for this batch.
    """
    degree, braid, start, end, progress, braid_state = args
    n = end - start
    for i, idx in enumerate(range(start, end)):
        if progress and n > 0 and (i % max(1, n // 10) == 0):
            pct = round(100.0 * i / n, 1)
            print(f"[worker {os.getpid()}] {pct}% of batch", end="\r")
        hit = check_assignment_for_braid(idx, degree, braid, braid_state)
        if hit is not None:
            return hit
    return None


# 7) parallel_try_sign_assignments
import multiprocessing as mp

def parallel_try_sign_assignments(
    degree,
    braid,
    braid_state,
    max_workers=None,
    chunk_size=1 << 14,   # tune: 16K indices per task by default
    include_flip=True,
    max_shifts=None,
    verbose=False,
):
    """
    Parallel version of your search routine.
    Strategy:
      - Iterate braid variants (original, flipped, shifts).
      - For each variant, split the index space into ranges and parallelize.
      - Return first (braid_variant, strand_signs) that validates, else False.
    """
    if max_workers is None:
        max_workers = max(1, (os.cpu_count() or 2) - 0)

    # Guard against pathological recursion in original version:
    original_tuple = tuple(braid)

    for variant_braid, meta in generate_braid_variants(braid, include_flip=include_flip, max_shifts=max_shifts):
        if tuple(variant_braid) == original_tuple and meta.get("shift", 0) == 0 and meta.get("flipped", False) is False:
            if verbose:
                print("Trying original braid...")
        elif meta.get("flipped", False):
            if verbose:
                print("Trying flipped braid...")
        else:
            if verbose:
                print(f"Trying cyclic shift by {meta['shift']}...")

        total = total_assignments_for_braid(braid_state)
        if verbose:
            print(f"Search space: {total:,} assignments")

        # Create tasks
        tasks = [
            (degree, variant_braid, start, end, verbose, braid_state)
            for (start, end) in split_ranges(total, chunk_size)
        ]

        # Run pool
        with mp.Pool(processes=max_workers) as pool:
            async_results = [pool.apply_async(worker_on_range, (t,)) for t in tasks]

            # gather until first success
            for r in async_results:
                out = r.get()
                if out is not None:
                    # found a solution; terminate the pool early
                    try:
                        pool.terminate()
                    finally:
                        pool.join()
                    return out

            # No success for this variant; ensure pool is closed
            pool.close()
            pool.join()

    # All variants exhausted
    return False


# 8) get_sign_assignment_parallel (drop-in wrapper with same schema as your original)
def get_sign_assignment_parallel(braid, degree=10, verbose=False, **kwargs):
    """
    Parallel version mirroring get_sign_assignment's output.
    kwargs are passed to parallel_try_sign_assignments (e.g., max_workers, chunk_size).
    """
    braid_type_ = braid_type(braid)
    braid_states = BraidStates(braid)

    if braid_type_ == BRAID_TYPE_HOMEGENOUS:
        # homogeneous: just return the initial strand_signs (same as your original)
        sign_assignment = braid_states.strand_signs
        return {
            'link_type': 'homogeneous',
            'inversion_data': sign_assignment,
            'braid': braid_states.braid,
            'degree': degree
        }

    elif braid_type_ == BRAID_TYPE_FIBERED:
        sol = parallel_try_sign_assignments(degree, braid, braid_states, verbose=verbose, **kwargs)
        if sol is False:
            return {'inversion_data': 'failure'}

        new_braid, sign_assignment = sol
        return {
            'link_type': 'fibered',
            'inversion_data': sign_assignment,
            'braid': new_braid,
            'degree': degree
        }

if __name__ == "__main__":

    print('L9n10_0')
    import time
    braid = [-1, 2, 3, -2, -1, -2, -1, -2, -3, 2]
    st = time.time()
    assignment = get_sign_assignment_parallel(
        braid=braid,
        degree=10,
        verbose=True,
        max_workers=None,   # or set an explicit int
        chunk_size=1<<14,   # tune if you like
        include_flip=True,
        max_shifts=None     # or a small int to limit shifts
    )
    end = time.time()
    print(assignment, end="\n\n")
    print(end-st)
