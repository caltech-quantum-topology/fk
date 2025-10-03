#!/usr/bin/env python3
"""
Example usage of print_symbolic_relations function
"""

from fkcompute.braidstates_links import BraidStates
from fkcompute.relations_links import full_reduce, print_symbolic_relations

# Example: Fibered knot with inversion data
braid = [-1, -2, -2, -2, -1, 2, 2, 2]
degree = 5
inversion_data = {0: [-1, -1, -1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1]}

print("=== EXAMPLE: Symbolic Relations for Fibered Knot ===\n")

# Set up braid states
braid_states = BraidStates(braid)
braid_states.strand_signs = inversion_data
braid_states.compute_matrices()
braid_states.generate_position_assignments()
braid_states.compute_r_matrices()  # Need this for get_state_relations()

# Get and reduce relations
all_relations = braid_states.get_state_relations()
relations = full_reduce(all_relations)

# Print symbolic relations
print_symbolic_relations(degree, relations, braid_states)
