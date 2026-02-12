"""
SignedBraid class combining topology with sign assignment.

A SignedBraid extends BraidTopology with sign assignment data,
crossing matrices, R-matrices, and relation generation.
"""

import itertools
from typing import List, Dict, Any, Optional, Tuple

from .topology import BraidTopology
from .types import ZERO_STATE, NUNITY_STATE
from .word import is_homogeneous_braid


class SignedBraid:
    """
    A braid with sign assignment data.

    Holds a BraidTopology reference and adds sign-dependent state:
    strand_signs, matrices, r_matrices, sign_assignment.

    Parameters
    ----------
    topology
        BraidTopology instance providing the pure topological data.
    """

    def __init__(self, topology: BraidTopology):
        self.topology = topology
        self.strand_signs: Dict[int, List[int]] = {c: [] for c in range(topology.n_components)}

        if is_homogeneous_braid(topology.braid):
            self._init_homogeneous_signs()

    def _init_homogeneous_signs(self):
        """Initialize signs for homogeneous braids."""
        topo = self.topology
        strand_signs_ = {0: 1}
        for i in range(1, topo.max_strand + 1):
            if i in topo.braid:
                strand_signs_[i] = 1
            elif -i in topo.braid:
                strand_signs_[i] = -1
            else:
                raise Exception(f'expected one of +{i} or -{i} to appear in the braid')
        self.strand_signs = {i: [] for i in range(topo.n_components)}
        for component in list(self.strand_signs.keys()):
            for pair in topo.strand_endpoints[component]:
                self.strand_signs[component].append(strand_signs_[pair[0][0]])
        self.compute_matrices()
        self.generate_position_assignments()
        self.compute_r_matrices()

    def __getattr__(self, name):
        """Delegate attribute access to topology for backward compatibility."""
        # Avoid infinite recursion during init
        if name == 'topology':
            raise AttributeError(name)
        return getattr(self.topology, name)

    def compute_r_matrices(self) -> Optional[int]:
        """Compute R-matrices for each crossing based on sign assignments."""
        topo = self.topology
        self.r_matrices = []
        for index in range(topo.n_crossings):
            crossing_sign = topo.crossing_signs[index]
            if crossing_sign == 1:
                if self.matrices[index][0][1] * self.matrices[index][1][0] > 0:
                    self.r_matrices.append('R1')
                elif self.matrices[index][0][1] == 1 and self.matrices[index][1][0] == -1:
                    self.r_matrices.append('R2')
                else:
                    return None
            elif crossing_sign == -1:
                if self.matrices[index][0][0] * self.matrices[index][1][1] > 0:
                    self.r_matrices.append('R4')
                elif self.matrices[index][0][0] == 1 and self.matrices[index][1][1] == -1:
                    self.r_matrices.append('R3')
                else:
                    return None
        return 0

    def validate(self) -> bool:
        """Validate the current sign assignment."""
        for index in range(self.topology.n_crossings):
            if self.matrices[index][0].count(1) != self.matrices[index][1].count(1):
                return False
        if self.compute_r_matrices() is not None:
            return True
        return False

    def none_count(self, matrix: List[List[Any]]) -> int:
        """Count the number of None values in a matrix."""
        return list(itertools.chain.from_iterable(matrix)).count(None)

    def sign_at_position(self, position: Tuple[int, int]) -> Optional[int]:
        """Get the sign at a given position in the crossing matrices."""
        topo = self.topology
        for index in range(topo.n_crossings):
            for index1 in range(2):
                for index2 in range(2):
                    if topo.loc[index][index1][index2] == position:
                        return self.matrices[index][index1][index2]
        return None

    def generate_strand_signs(self):
        """Generate strand signs from the current assignment."""
        self.strand_signs = []
        for index in range(1, self.topology.n_s):
            self.strand_signs.append(self.sign_assignment[self.topology.strand_endpoints[index][0]])

    def generate_position_assignments(self):
        """Generate position assignments from strand signs."""
        topo = self.topology
        self.sign_assignment = dict()
        for component in range(topo.n_components):
            for index in range(topo.n_s[component]):
                sign = self.sign_at_position(topo.strand_locations[component][index][0])
                for index_ in range(len(topo.strand_locations[component][index])):
                    self.sign_assignment[topo.strand_locations[component][index][index_]] = sign

    def compute_matrices(self):
        """Compute crossing matrices from strand signs."""
        topo = self.topology
        self.matrices = [[[None, None], [None, None]] for _ in range(topo.n_crossings)]
        for component in range(topo.n_components):
            for index in range(topo.n_s[component]):
                if topo.strand_types[component][index][0] == 'L':
                    bottom = 0
                elif topo.strand_types[component][index][0] == 'R':
                    bottom = 1
                if topo.strand_types[component][index][1] == 'L':
                    top = 0
                elif topo.strand_types[component][index][1] == 'R':
                    top = 1
                self.matrices[topo.endpoint_crossing_indices[component][index][0]][0][bottom] = self.strand_signs[component][index]
                self.matrices[topo.endpoint_crossing_indices[component][index][1]][1][top] = self.strand_signs[component][index]

    def update(self, position: Tuple[int, int], value: int):
        """Update a position in the crossing matrices."""
        topo = self.topology
        for index in range(topo.n_s):
            for index_ in range(2):
                if topo.strand_endpoints[index][index_] == position:
                    position = topo.strand_endpoints[index][(index_ + 1) % 2]
                    break
        for index in range(topo.n_crossings):
            for index1 in range(2):
                for index2 in range(2):
                    if topo.loc[index][index1][index2] == position:
                        if self.matrices[index][index1][index2] == -value:
                            raise Exception("Tried to overwrite R-matrix sign information!")
                        else:
                            self.matrices[index][index1][index2] = value

    def propagate(self):
        """Propagate sign constraints through the crossing matrices."""
        topo = self.topology
        index = 0
        while index < topo.n_crossings:
            reset = False
            count = self.none_count(self.matrices[index])
            if count == 2 or count == 1:
                for index1 in range(2):
                    for index2 in range(2):
                        if self.matrices[index][0][index1] == 1 and self.matrices[index][1][index2] == -1:
                            self.matrices[index][0][(index1 + 1) % 2] = -1
                            self.matrices[index][1][(index2 + 1) % 2] = 1
                            position1 = topo.loc[index][0][(index1 + 1) % 2]
                            self.update(position1, -1)
                            position2 = topo.loc[index][1][(index2 + 1) % 2]
                            self.update(position2, 1)
                            reset = True
                        elif self.matrices[index][0][index1] == -1 and self.matrices[index][1][index2] == 1:
                            self.matrices[index][0][(index1 + 1) % 2] = 1
                            self.matrices[index][1][(index2 + 1) % 2] = -1
                            position1 = topo.loc[index][0][(index1 + 1) % 2]
                            self.update(position1, 1)
                            position2 = topo.loc[index][1][(index2 + 1) % 2]
                            self.update(position2, -1)
                            reset = True
            elif count == 1:
                for index1 in range(2):
                    if self.matrices[index1][0] == [1, 1]:
                        for index2 in range(2):
                            if self.matrices[(index1 + 1) % 2][1][index2] == 1:
                                self.matrices[(index1 + 1) % 2][1][(index2 + 1) % 2] = 1
                                position = topo.loc[(index1 + 1) % 2][1][(index2 + 1) % 2]
                                self.update(position, 1)
                                reset = True
                    elif self.matrices[index1][0] == [-1, -1]:
                        for index2 in range(2):
                            if self.matrices[(index1 + 1) % 2][1][index2] == -1:
                                self.matrices[(index1 + 1) % 2][1][(index2 + 1) % 2] = -1
                                position = topo.loc[(index1 + 1) % 2][1][(index2 + 1) % 2]
                                self.update(position, -1)
                                reset = True
            if reset:
                index = 0
            else:
                index += 1

    def get_state_relations(self) -> List:
        """Generate state relations for the current braid and sign assignment."""
        from ..constraints.relations import Zero, Nunity, Alias, Conservation, Leq

        topo = self.topology
        relations = []
        if self.strand_signs[0][0] == 1:
            relations.append(Zero(topo.get_state((0, 0))))
            relations.append(Zero(topo.get_state((0, len(topo.braid)))))
        else:
            relations.append(Nunity(topo.get_state((0, 0))))
            relations.append(Nunity(topo.get_state((0, len(topo.braid)))))

        for i in topo.strands:
            relations.append(Alias(topo.get_state((i, len(topo.braid))), topo.get_state((i, 0))))

        for (i, j) in topo.state_info.keys():
            if self.sign_assignment[(i, j)] == 1:
                relations.append(Leq(ZERO_STATE, (i, j)))
            elif self.sign_assignment[(i, j)] == -1:
                relations.append(Leq((i, j), NUNITY_STATE))
            else:
                raise ValueError(f"The sign at position {(i, j)} is not 1 or -1!")

        for j, gen in enumerate(topo.braid):
            if self.r_matrices[j] == 'R1':
                relations.append(Leq(topo.get_state((abs(gen), j + 1)), topo.get_state((abs(gen) - 1, j))))
                relations.append(Leq(topo.get_state((abs(gen), j)), topo.get_state((abs(gen) - 1, j + 1))))
            elif self.r_matrices[j] == 'R4':
                relations.append(Leq(topo.get_state((abs(gen) - 1, j)), topo.get_state((abs(gen), j + 1))))
                relations.append(Leq(topo.get_state((abs(gen) - 1, j + 1)), topo.get_state((abs(gen), j))))

            relations.append(Conservation(
                [topo.get_state((abs(gen) - 1, j)), topo.get_state((abs(gen), j))],
                [topo.get_state((abs(gen) - 1, j + 1)), topo.get_state((abs(gen), j + 1))]
            ))
        return relations

    def reduced_relations(self) -> List:
        """Get the reduced relations for this braid."""
        from ..constraints.reduction import full_reduce
        return full_reduce(self.get_state_relations())

    def free_variables(self) -> List:
        """Get the free variables after reduction."""
        from ..constraints.reduction import free_variables
        return [x for x in free_variables(self.reduced_relations()) if x != ZERO_STATE and x != NUNITY_STATE]

    def load_sign_data(self, sign_data: List[int]):
        """Load sign data into the braid states."""
        from ..constraints.reduction import full_reduce

        topo = self.topology
        for component in range(topo.n_components):
            self.strand_signs[component] = sign_data[:topo.n_s[component]]
            sign_data = sign_data[topo.n_s[component]:]
        self.compute_matrices()
        if self.validate():
            self.generate_position_assignments()
            all_relations = self.get_state_relations()
            relations = full_reduce(all_relations)
            return relations
        return False
