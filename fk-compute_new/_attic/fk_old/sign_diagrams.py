from braidstates import BraidStates
from relations import full_reduce
from braid_ilp import check_sign_assignment
from classification import is_homogeneous_braid

BRAID_TYPE_HOMEGENOUS = 0
BRAID_TYPE_FIBERED = 1

def braid_type(braid):
    if is_homogeneous_braid(braid):
        return BRAID_TYPE_HOMEGENOUS
    else:
        return BRAID_TYPE_FIBERED

def load_sign_data(
        bs : BraidStates,
        sign_data=None,
        compute_r_matrices=True
    ) -> bool:
    if is_homogeneous_braid(bs.braid):
            bs.strand_signs_ = {0:1}
            for i in range(1,bs.max_strand+1):
                if i in bs.braid:
                    bs.strand_signs_[i] = 1
                elif -i in bs.braid:
                    bs.strand_signs_[i] = -1
                else:
                    raise Exception('expected one of +{i} or -{i} to appear in the braid')
            bs.strand_signs = {i:[] for i in range(bs.n_components)}
            for component in list(bs.strand_signs.keys()):
                for pair in bs.strand_endpoints[component]:
                    bs.strand_signs[component].append(bs.strand_signs_[pair[0][0]])
            bs.compute_matrices()
            bs.generate_position_assignments()
            if compute_r_matrices:
                bs.compute_r_matrices()
            return True
    if isinstance(sign_data, dict):
        values = sign_data.values()
        sign_data = []
        for v in values:
            sign_data += v
    if sign_data is None or len(sign_data) != bs.n_s_total:
        raise Exception('invalid sign data for non-homogeneous braid')
    for component in range(bs.n_components):
        bs.strand_signs[component] = sign_data[:bs.n_s[component]]
        sign_data = sign_data[bs.n_s[component]:]
    bs.compute_matrices()
    if bs.validate():
        bs.generate_position_assignments()
        if compute_r_matrices:
            bs.compute_r_matrices()
        return True
    return False

def flip(braid):
    n = BraidStates(braid).n_strands
    return list(map(lambda x: (2 * (x > 0) - 1) * (n - abs(x)), braid))

def try_sign_assignments(degree, braid, weight=-1, max_depth=-1, flipped=False, verbose=False, original_braid = []):

    braid_states = BraidStates(braid)

    # look into making the sign diagram search more modular, isolating the ILP and bound checks, while making the individual knot theory modules more clean and efficient

    # also, and pressingly, identify why the code explodes your computer sometimes

    if braid == original_braid:
        return False
    if original_braid == []:
        original_braid = braid[:]

    index = 0
    n = 2**braid_states.n_s_total
    for index in range(n):
        print(f"finding inversion data... {round(100 * index/n, 1)}%", end='\r')
        sign_assignment = list(map(lambda x: 2 * int(x) - 1, list(str(bin(index))[2:].zfill(braid_states.n_s_total))))
        if load_sign_data(braid_states, sign_assignment, compute_r_matrices=False):
            all_relations = braid_states.get_state_relations()
            relations = full_reduce(all_relations)
            out = check_sign_assignment(degree, relations, braid_states, weight)
            if out is not None:
                return braid, braid_states.strand_signs

    # if we're not successful on the first pass, we need to try flipping/shifting the braid (or some other operations)

    if not flipped:
        return try_sign_assignments(degree=degree, braid=flip(braid), max_depth=max_depth, flipped=True, original_braid=original_braid)
    else:
        braid = flip(braid)
        return try_sign_assignments(degree=degree, braid=braid[1:] + [braid[0]], max_depth=max_depth, flipped=False, original_braid=original_braid)

def get_sign_assignment(braid, degree=10, verbose=False):
    braid_type_ = braid_type(braid)
    braid_states = BraidStates(braid)
    if braid_type_ == BRAID_TYPE_HOMEGENOUS:
        sign_assignment = braid_states.strand_signs
        return {
            'link_type' : 'homogeneous',
            'inversion_data': sign_assignment,
            'braid': braid_states.braid,
            'degree': degree
        }
    elif braid_type_ == BRAID_TYPE_FIBERED:
        sol = try_sign_assignments(degree, braid, verbose=verbose)
        if sol == False:
            return {
                'inversion_data': 'failure'
            }
        new_braid, sign_assignment = sol
        return {
            'link_type' : 'fibered',
            'inversion_data': sign_assignment,
            'braid': new_braid,
            'degree': degree
        }
