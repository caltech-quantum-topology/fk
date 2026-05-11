"""Compare regenerated ILP CSV inequalities to the existing Unknown_161658_ilp.csv."""

import sys
sys.path.insert(0, '/home/davep96/work/projects/fk_computation/fk/fk-compute/src')

import numpy as np
from fkcompute import BraidStates
from fkcompute.domain.constraints.reduction import full_reduce
from fkcompute.solver.ilp import _czech_sign_assignment, ilp
from fkcompute.domain.constraints.relations import _sort_any

braid = [2, -3, 2, 3, 3, 3, 3, -2, -1, 2, -3, -1]
inversion = {
    0: [-1, -1, 1, 1, 1, 1, 1, -1],
    1: [-1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, 1, -1, -1, -1],
}
degree = 3

braid_states = BraidStates(braid)
braid_states.strand_signs = inversion
braid_states.compute_matrices()
braid_states.compute_r_matrices()
braid_states.generate_position_assignments()

all_relations = braid_states.get_state_relations()
relations = full_reduce(all_relations)

check = _czech_sign_assignment(degree, relations, braid_states)
if check is None:
    print("ERROR: No valid sign assignment found")
    sys.exit(1)

# Regenerate the CSV string
generated = ilp(degree, relations, braid_states)
generated_lines = generated.strip().split('\n')

# Read the existing CSV
with open('Unknown_161658_ilp.csv') as f:
    existing_lines = [l.rstrip() for l in f.readlines()]

# Strip trailing slashes/whitespace for comparison
def parse_tableau(lines):
    rows = []
    for l in lines:
        l = l.rstrip(',').strip()
        if l and l != '/':
            rows.append([float(x) for x in l.split(',')])
    return rows

# Split both into sections (header, criteria, inequalities, assignments)
def split_sections(lines):
    sections = []
    current = []
    for l in lines:
        l = l.rstrip()
        if l == '/':
            sections.append(current)
            current = []
        else:
            current.append(l)
    sections.append(current)
    return sections

gen_sections = split_sections(generated_lines)
ext_sections = split_sections(existing_lines)

print(f"Generated sections: {len(gen_sections)}, Existing sections: {len(ext_sections)}")

# Section 0: header (degree, components, writhe, braid, closed strands, crossing components)
# Section 1: criteria
# Section 2: inequalities
# Section 3: assignments

section_names = ['header', 'criteria', 'inequalities', 'assignments']

for si in range(min(len(gen_sections), len(ext_sections))):
    name = section_names[si] if si < len(section_names) else f'section_{si}'
    gen_rows = [l.rstrip(',') for l in gen_sections[si] if l.strip()]
    ext_rows = [l.rstrip(',') for l in ext_sections[si] if l.strip()]

    if si == 0:
        # Just compare header lines
        if gen_rows != ext_rows:
            print(f"\n[{name}] MISMATCH:")
            print(f"  Generated: {gen_rows}")
            print(f"  Existing:  {ext_rows}")
        else:
            print(f"[{name}] OK")
        continue

    # Parse as numeric
    try:
        gen_mat = np.array([[float(x) for x in r.split(',')] for r in gen_rows]) if gen_rows else np.zeros((0,0))
        ext_mat = np.array([[float(x) for x in r.split(',')] for r in ext_rows]) if ext_rows else np.zeros((0,0))
    except Exception as e:
        print(f"[{name}] Parse error: {e}")
        continue

    print(f"\n[{name}] Generated shape: {gen_mat.shape}, Existing shape: {ext_mat.shape}")

    if gen_mat.shape != ext_mat.shape:
        print(f"  SHAPE MISMATCH")
        if gen_mat.size > 0:
            print(f"  Generated rows:")
            for row in gen_mat:
                print(f"    {list(row)}")
        if ext_mat.size > 0:
            print(f"  Existing rows:")
            for row in ext_mat:
                print(f"    {list(row)}")
        continue

    # Check if rows match (possibly reordered)
    gen_set = set(tuple(r) for r in gen_mat.tolist())
    ext_set = set(tuple(r) for r in ext_mat.tolist())

    missing_from_existing = gen_set - ext_set
    missing_from_generated = ext_set - gen_set

    if not missing_from_existing and not missing_from_generated:
        if np.allclose(gen_mat, ext_mat):
            print(f"  IDENTICAL (same order)")
        else:
            print(f"  MATCH (different order)")
    else:
        print(f"  MISMATCH")
        if missing_from_existing:
            print(f"  In generated but NOT in existing ({len(missing_from_existing)} rows):")
            for r in sorted(missing_from_existing):
                print(f"    {list(r)}")
        if missing_from_generated:
            print(f"  In existing but NOT in generated ({len(missing_from_generated)} rows):")
            for r in sorted(missing_from_generated):
                print(f"    {list(r)}")
