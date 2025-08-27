from braidstates_links import BraidStates
from relations_links import full_reduce, check_sign_assignment
from braids import is_homogeneous_braid
import json

# --- classification constants -------------------------------------------------
BRAID_TYPE_HOMEGENOUS = 0
BRAID_TYPE_FIBERED = 1

def braid_type(braid):
    return BRAID_TYPE_HOMEGENOUS if is_homogeneous_braid(braid) else BRAID_TYPE_FIBERED

# --- simple transforms --------------------------------------------------------
def flip(braid):
    """Mirror/flip the braid indices using the number of strands given by BraidStates."""
    n = BraidStates(braid).n_strands
    return [(2 * (x > 0) - 1) * (n - abs(x)) for x in braid]

def rotate(braid):
    """Rotate braid left by one position."""
    return braid[1:] + [braid[0]]

def generate_braid_variants(braid, include_flip=True, max_shifts=None):
    """
    Yield braid variants in order:
      shift=0: braid, braid flipped (if include_flip=True)
      shift=1: rotated once, rotated once flipped
      shift=2: rotated twice, rotated twice flipped
      ...
    Stops after len(braid) shifts, or earlier if max_shifts is set.
    Yields: (variant_braid, {"flipped": bool, "shift": int})
    """
    limit = len(braid) if max_shifts is None else min(len(braid), max_shifts)
    current = braid[:]
    for shift in range(limit):
        yield current, {"flipped": False, "shift": shift}
        if include_flip:
            yield flip(current), {"flipped": True, "shift": shift}
        current = rotate(current)

# --- sign assignment utilities -----------------------------------------------
def sign_assignment_from_index(index, braid_states):
    """
    Convert an integer index to per-component sign lists.
    Uses binary expansion of 'index' over n_s_total bits and maps bit -> sign via 2*bit-1.
    Returns: list[list[int]] with shape [n_components][n_s[component]]
    """
    bits = list(map(int, list(bin(index)[2:].zfill(braid_states.n_s_total))))
    signs = [2 * b - 1 for b in bits]  # 0 -> -1, 1 -> +1

    out, cursor = [], 0
    for component in range(braid_states.n_components):
        k = braid_states.n_s[component]
        out.append(signs[cursor:cursor + k])
        cursor += k
    return out

def check_assignment_for_braid(index, degree, braid, braid_state):
    """
    Try a single sign-assignment (by index) for 'braid'.
    If valid, return (braid, strand_signs); else return None.
    """

    # fill strand signs from index
    per_component = sign_assignment_from_index(index, braid_state)
    for component in range(braid_state.n_components):
        braid_state.strand_signs[component] = per_component[component]

    # compute & validate
    braid_state.compute_matrices()
    if not braid_state.validate():
        return None

    braid_state.generate_position_assignments()
    all_relations = braid_state.get_state_relations()
    relations = full_reduce(all_relations)
    out = check_sign_assignment(degree, relations, braid_state)

    if out is not None:
        return braid, braid_state.strand_signs
    return None

# --- serial search ------------------------------------------------------------
def try_sign_assignments_serial(
    degree,
    braid,
    braid_state,
    include_flip=True,
    max_shifts=None,
    verbose=False,
):
    """
    Serial (non-parallel) search for a consistent sign assignment.
    Iterates braid variants in the specified order, then enumerates all 2^{n_s_total}
    sign assignments for each variant until one validates.
    Returns: (braid_variant, strand_signs) on success, or False on failure.
    """
    variants_total = 2*len(braid)
    variant_idx = -1
    for variant_braid, meta in generate_braid_variants(braid, include_flip=include_flip, max_shifts=max_shifts):
        variant_idx +=1
        if verbose:
            head = "flipped" if meta["flipped"] else "original"
            if meta["shift"] == 0:
                print(f"Trying {head} braid (shift={meta['shift']})...")
            else:
                print(f"Trying rotation shift={meta['shift']} ({head})...")

        
        total = 2 ** braid_state.n_s_total
        if verbose:
            print(f"  Search space: {total:,} assignments")

        # serial enumeration
        for idx in range(total):
            if verbose:
                pct = round(100.0 * (variant_idx/variants_total + idx / (variants_total*total)), 2)
                print(f"    {pct}%", end="\r")

            hit = check_assignment_for_braid(idx, degree, variant_braid, braid_state)
            if hit is not None:
                if verbose:
                    print("\n  ✓ Found valid assignment.")
                return hit

        if verbose:
            print("  No valid assignment for this variant.\n")

    # No success across all variants
    return False

def get_sign_assignment_serial(braid, degree=10, verbose=False, **kwargs):
    """
    Non-parallel version mirroring your original 'get_sign_assignment' return schema.
    kwargs are accepted for symmetry but ignored here (kept for drop-in replacement flexibility).
    """
    t = braid_type(braid)
    bs = BraidStates(braid)

    if t == BRAID_TYPE_HOMEGENOUS:
        # homogeneous: just return the initial strand_signs
        sign_assignment = bs.strand_signs
        return {
            'link_type': 'homogeneous',
            'inversion_data': sign_assignment,
            'braid': bs.braid,
            'degree': degree
        }

    # fibered / non-homogeneous: run serial search
    sol = try_sign_assignments_serial(degree, braid, bs, verbose=verbose)
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
    import time
    braid = [-1, 2, 3, -2, -1, -2, -1, -2, -3, 2]
    st = time.time()
    assignment = get_sign_assignment_serial(
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

