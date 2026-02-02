"""
BraidStates class for braid state management.

The main class BraidStates represents the combined data of:
- a braid
- the states of the braid
- known relations among states
"""

import string
import itertools
import copy
from typing import List, Dict, Any, Optional, Tuple

from .types import ZERO_STATE, NUNITY_STATE
from .word import is_homogeneous_braid


class BraidStates:
    """
    Represents a braid and its associated state information.

    This class tracks the topology of a braid closure, including:
    - State locations on the braid diagram
    - Component structure (which strands form which components)
    - Sign assignments at crossings
    - Relations between states

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.
    """

    def __init__(self, braid: List[int]):
        self.braid = braid
        self.max_strand = max(abs(x) for x in braid)
        self.strands = list(range(0, self.max_strand + 1))
        self.n_strands = self.max_strand + 1
        self.n_crossings = len(braid)
        self.braid_group_generators = list(range(1, self.max_strand + 1))
        self.crossing_signs = []
        for g in braid:
            if g < 0:
                self.crossing_signs.append(-1)
            else:
                self.crossing_signs.append(+1)
        self.writhe = sum(self.crossing_signs)

        self.state_locations = [(i, j) for i in self.strands for j in range(0, self.n_crossings + 1)]
        self.top_input_state_locations = [(abs(braid[index]) - 1, index) for index in range(self.n_crossings)]
        self.bottom_input_state_locations = [(x + 1, y) for (x, y) in self.top_input_state_locations]
        self.top_output_state_locations = [(x, y + 1) for (x, y) in self.top_input_state_locations]
        self.bottom_output_state_locations = [(x, y + 1) for (x, y) in self.bottom_input_state_locations]
        self.loc = [[[(x, y + 1), (x + 1, y + 1)], [(x, y), (x + 1, y)]] for (x, y) in self.top_input_state_locations]

        state_info = {}

        state_at_location = {}

        for strand in self.strands:
            state = (strand, 0)
            state_info[state] = {}
            state_at_location[(strand, 0)] = (strand, 0)
            for j in range(1, self.n_crossings + 1):
                if abs(braid[j - 1]) == strand or abs(braid[j - 1]) == strand + 1:
                    state = (strand, j)
                    state_info[state] = {}
                state_at_location[(strand, j)] = state

        for i, (state, info) in enumerate(state_info.items()):
            info['label'] = str(i) if len(state_info) > 26 else string.ascii_lowercase[i]

        self.state_info = state_info
        self._state_at_location = state_at_location
        self.states = list(sorted(set(state_at_location.values())))

        self.state_equivalence_classes = {}
        for state in self.states:
            self.state_equivalence_classes[state] = []
        for key in state_at_location.keys():
            self.state_equivalence_classes[state_at_location[key]].append(key)

        self.component_locations = []
        self.strand_locations = dict()
        index = 0
        component = 0
        while index < self.n_strands:
            current_loc = [index, 0]
            if tuple(current_loc) not in list(itertools.chain.from_iterable(self.component_locations)):
                done = False
                prefix = []
                self.strand_locations[component] = []
                current_list = []
                current_list_strand = []
                while not done:
                    current_list.append(tuple(current_loc))
                    current_list_strand.append(tuple(current_loc))
                    if tuple(current_loc) in self.bottom_input_state_locations:
                        self.strand_locations[component].append(prefix + current_list_strand + ['R'])
                        current_loc = [current_loc[0] - 1, current_loc[1] + 1]
                        current_list_strand = []
                        prefix = ['L']
                    elif tuple(current_loc) in self.top_input_state_locations:
                        self.strand_locations[component].append(prefix + current_list_strand + ['L'])
                        current_loc = [current_loc[0] + 1, current_loc[1] + 1]
                        current_list_strand = []
                        prefix = ['R']
                    elif current_loc[1] == self.n_crossings:
                        current_loc = [current_loc[0], 0]
                    else:
                        current_loc = [current_loc[0], current_loc[1] + 1]
                    if tuple(current_loc) in current_list:
                        self.strand_locations[component][0] = prefix + current_list_strand + self.strand_locations[component][0]
                        done = True
                if current_list != []:
                    self.component_locations.append(current_list)
                    component += 1
            index += 1

        self.component_locations_dict = {index: list_ for (index, list_) in enumerate(self.component_locations)}
        self.n_components = len(self.component_locations)
        self.top_crossing_components = [self.get_component(location) for location in self.top_input_state_locations]
        self.bottom_crossing_components = [self.get_component(location) for location in self.bottom_input_state_locations]
        self.closed_strand_components = [self.get_component(location) for location in [(index, 0) for index in range(self.n_strands)]]

        self.strand_types = {key: [[x[0], x[-1]] for x in value] for (key, value) in self.strand_locations.items()}
        self.strand_locations = {key: [x[1:-1] for x in value] for (key, value) in self.strand_locations.items()}
        self.strand_endpoints = {key: [[x[0], x[-1]] for x in value] for (key, value) in self.strand_locations.items()}
        self.endpoint_crossing_indices = {key: [[x[0][1] - 1, x[-1][1]] for x in value] for (key, value) in self.strand_endpoints.items()}
        self.n_s = {key: len(value) for (key, value) in self.strand_locations.items()}
        self.n_s_total = sum(list(self.n_s.values()))
        self.strand_signs: Dict[int, List[int]] = {c: [] for c in range(self.n_components)}

        if is_homogeneous_braid(braid):
            self.strand_signs_ = {0: 1}
            for i in range(1, self.max_strand + 1):
                if i in braid:
                    self.strand_signs_[i] = 1
                elif -i in braid:
                    self.strand_signs_[i] = -1
                else:
                    raise Exception(f'expected one of +{i} or -{i} to appear in the braid')
            self.strand_signs = {i: [] for i in range(self.n_components)}
            for component in list(self.strand_signs.keys()):
                for pair in self.strand_endpoints[component]:
                    self.strand_signs[component].append(self.strand_signs_[pair[0][0]])
            self.compute_matrices()
            self.generate_position_assignments()
            self.compute_r_matrices()

    def compute_r_matrices(self) -> Optional[int]:
        """Compute R-matrices for each crossing based on sign assignments."""
        self.r_matrices = []
        for index in range(self.n_crossings):
            crossing_sign = self.crossing_signs[index]
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
        for index in range(self.n_crossings):
            if self.matrices[index][0].count(1) != self.matrices[index][1].count(1):
                return False
        if self.compute_r_matrices() is not None:
            return True
        return False

    def get_component(self, location: Tuple[int, int]) -> Optional[int]:
        """Get the component index for a given location."""
        for (key, value) in self.component_locations_dict.items():
            if location in value:
                return key
        return None

    def none_count(self, matrix: List[List[Any]]) -> int:
        """Count the number of None values in a matrix."""
        return list(itertools.chain.from_iterable(matrix)).count(None)

    def sign_at_position(self, position: Tuple[int, int]) -> Optional[int]:
        """Get the sign at a given position in the crossing matrices."""
        for index in range(self.n_crossings):
            for index1 in range(2):
                for index2 in range(2):
                    if self.loc[index][index1][index2] == position:
                        return self.matrices[index][index1][index2]
        return None

    def generate_strand_signs(self):
        """Generate strand signs from the current assignment."""
        self.strand_signs = []
        for index in range(1, self.n_s):
            self.strand_signs.append(self.sign_assignment[self.strand_endpoints[index][0]])

    def generate_position_assignments(self):
        """Generate position assignments from strand signs."""
        self.sign_assignment = dict()
        for component in range(self.n_components):
            for index in range(self.n_s[component]):
                sign = self.sign_at_position(self.strand_locations[component][index][0])
                for index_ in range(len(self.strand_locations[component][index])):
                    self.sign_assignment[self.strand_locations[component][index][index_]] = sign

    def compute_matrices(self):
        """Compute crossing matrices from strand signs."""
        self.matrices = [[[None, None], [None, None]] for _ in range(self.n_crossings)]
        for component in range(self.n_components):
            for index in range(self.n_s[component]):
                if self.strand_types[component][index][0] == 'L':
                    bottom = 0
                elif self.strand_types[component][index][0] == 'R':
                    bottom = 1
                if self.strand_types[component][index][1] == 'L':
                    top = 0
                elif self.strand_types[component][index][1] == 'R':
                    top = 1
                self.matrices[self.endpoint_crossing_indices[component][index][0]][0][bottom] = self.strand_signs[component][index]
                self.matrices[self.endpoint_crossing_indices[component][index][1]][1][top] = self.strand_signs[component][index]

    def update(self, position: Tuple[int, int], value: int):
        """Update a position in the crossing matrices."""
        for index in range(self.n_s):
            for index_ in range(2):
                if self.strand_endpoints[index][index_] == position:
                    position = self.strand_endpoints[index][(index_ + 1) % 2]
                    break
        for index in range(self.n_crossings):
            for index1 in range(2):
                for index2 in range(2):
                    if self.loc[index][index1][index2] == position:
                        if self.matrices[index][index1][index2] == -value:
                            raise Exception("Tried to overwrite R-matrix sign information!")
                        else:
                            self.matrices[index][index1][index2] = value

    def propagate(self):
        """Propagate sign constraints through the crossing matrices."""
        index = 0
        while index < self.n_crossings:
            reset = False
            count = self.none_count(self.matrices[index])
            if count == 2 or count == 1:
                for index1 in range(2):
                    for index2 in range(2):
                        if self.matrices[index][0][index1] == 1 and self.matrices[index][1][index2] == -1:
                            self.matrices[index][0][(index1 + 1) % 2] = -1
                            self.matrices[index][1][(index2 + 1) % 2] = 1
                            position1 = self.loc[index][0][(index1 + 1) % 2]
                            self.update(position1, -1)
                            position2 = self.loc[index][1][(index2 + 1) % 2]
                            self.update(position2, 1)
                            reset = True
                        elif self.matrices[index][0][index1] == -1 and self.matrices[index][1][index2] == 1:
                            self.matrices[index][0][(index1 + 1) % 2] = 1
                            self.matrices[index][1][(index2 + 1) % 2] = -1
                            position1 = self.loc[index][0][(index1 + 1) % 2]
                            self.update(position1, 1)
                            position2 = self.loc[index][1][(index2 + 1) % 2]
                            self.update(position2, -1)
                            reset = True
            elif count == 1:
                for index1 in range(2):
                    if self.matrices[index1][0] == [1, 1]:
                        for index2 in range(2):
                            if self.matrices[(index1 + 1) % 2][1][index2] == 1:
                                self.matrices[(index1 + 1) % 2][1][(index2 + 1) % 2] = 1
                                position = self.loc[(index1 + 1) % 2][1][(index2 + 1) % 2]
                                self.update(position, 1)
                                reset = True
                    elif self.matrices[index1][0] == [-1, -1]:
                        for index2 in range(2):
                            if self.matrices[(index1 + 1) % 2][1][index2] == -1:
                                self.matrices[(index1 + 1) % 2][1][(index2 + 1) % 2] = -1
                                position = self.loc[(index1 + 1) % 2][1][(index2 + 1) % 2]
                                self.update(position, -1)
                                reset = True
            if reset:
                index = 0
            else:
                index += 1

    def get_state(self, location: Tuple[int, int]) -> Optional[Tuple[int, int]]:
        """Get the state at a given location."""
        return self._state_at_location.get(location)

    def label_from_location(self, location) -> str:
        """Get the label for a location."""
        if type(location) != tuple:
            return str(location)
        return self.state_info[self.get_state(location)].get('label') or str(location)

    def get_state_relations(self) -> List:
        """Generate state relations for the current braid and sign assignment."""
        # Import here to avoid circular imports
        from ..constraints.relations import Zero, Nunity, Alias, Conservation, Leq

        relations = []
        if self.strand_signs[0][0] == 1:
            relations.append(Zero(self.get_state((0, 0))))
            relations.append(Zero(self.get_state((0, len(self.braid)))))
        else:
            relations.append(Nunity(self.get_state((0, 0))))
            relations.append(Nunity(self.get_state((0, len(self.braid)))))

        for i in self.strands:
            relations.append(Alias(self.get_state((i, len(self.braid))), self.get_state((i, 0))))

        for (i, j) in self.state_info.keys():
            if self.sign_assignment[(i, j)] == 1:
                relations.append(Leq(ZERO_STATE, (i, j)))
            elif self.sign_assignment[(i, j)] == -1:
                relations.append(Leq((i, j), NUNITY_STATE))
            else:
                raise ValueError(f"The sign at position {(i, j)} is not 1 or -1!")

        for j, gen in enumerate(self.braid):
            if self.r_matrices[j] == 'R1':
                relations.append(Leq(self.get_state((abs(gen), j + 1)), self.get_state((abs(gen) - 1, j))))
                relations.append(Leq(self.get_state((abs(gen), j)), self.get_state((abs(gen) - 1, j + 1))))
            elif self.r_matrices[j] == 'R4':
                relations.append(Leq(self.get_state((abs(gen) - 1, j)), self.get_state((abs(gen), j + 1))))
                relations.append(Leq(self.get_state((abs(gen) - 1, j + 1)), self.get_state((abs(gen), j))))

            relations.append(Conservation(
                [self.get_state((abs(gen) - 1, j)), self.get_state((abs(gen), j))],
                [self.get_state((abs(gen) - 1, j + 1)), self.get_state((abs(gen), j + 1))]
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

    def Lobb(self, L: int) -> bool:
        """Check the Lobb bound."""
        lower = sum([-x not in self.braid for x in range(1, self.n_strands)])
        upper = sum([x not in self.braid for x in range(1, self.n_strands)])
        U = self.writhe + self.n_strands + 1 - 2 * lower
        Delta = self.n_strands + 1 - upper - lower
        return L <= U and L >= U - 2 * Delta

    def a_resolution_components(self) -> Tuple[int, List]:
        """Compute the A-resolution components."""
        unvisited_locations = copy.deepcopy(self.state_locations)
        component_locations = []

        while unvisited_locations:
            current_loc = unvisited_locations[0]
            current_list = []
            up = True
            while tuple(current_loc) not in current_list:
                old_location = copy.deepcopy(current_loc)
                if up:
                    if tuple(current_loc) in self.top_input_state_locations:
                        if self.crossing_signs[self.top_input_state_locations.index(tuple(current_loc))] == 1:
                            current_loc = [current_loc[0], current_loc[1] + 1]
                        else:
                            current_loc = [current_loc[0] + 1, current_loc[1]]
                            up = False
                    elif tuple(current_loc) in self.bottom_input_state_locations:
                        if self.crossing_signs[self.bottom_input_state_locations.index(tuple(current_loc))] == 1:
                            current_loc = [current_loc[0], current_loc[1] + 1]
                        else:
                            current_loc = [current_loc[0] - 1, current_loc[1]]
                            up = False
                    elif current_loc[1] == self.n_crossings:
                        current_loc = [current_loc[0], 0]
                    else:
                        current_loc = [current_loc[0], current_loc[1] + 1]
                else:
                    if tuple(current_loc) in self.top_output_state_locations:
                        if self.crossing_signs[self.top_output_state_locations.index(tuple(current_loc))] == 1:
                            current_loc = [current_loc[0], current_loc[1] - 1]
                        else:
                            current_loc = [current_loc[0] + 1, current_loc[1]]
                            up = True
                    elif tuple(current_loc) in self.bottom_output_state_locations:
                        if self.crossing_signs[self.bottom_output_state_locations.index(tuple(current_loc))] == 1:
                            current_loc = [current_loc[0], current_loc[1] - 1]
                        else:
                            current_loc = [current_loc[0] - 1, current_loc[1]]
                            up = True
                    elif current_loc[1] == 0:
                        current_loc = [current_loc[0], self.n_crossings]
                    else:
                        current_loc = [current_loc[0], current_loc[1] - 1]
                current_list.append(tuple(old_location))
                unvisited_locations.remove(tuple(old_location))
            component_locations.append(current_list)

        return len(component_locations), component_locations

    def b_resolution_components(self) -> Tuple[int, List]:
        """Compute the B-resolution components."""
        unvisited_locations = copy.deepcopy(self.state_locations)
        component_locations = []

        while unvisited_locations:
            current_loc = unvisited_locations[0]
            current_list = []
            up = True
            while tuple(current_loc) not in current_list:
                old_location = copy.deepcopy(current_loc)
                if up:
                    if tuple(current_loc) in self.top_input_state_locations:
                        if self.crossing_signs[self.top_input_state_locations.index(tuple(current_loc))] == 1:
                            current_loc = [current_loc[0] + 1, current_loc[1]]
                            up = False
                        else:
                            current_loc = [current_loc[0], current_loc[1] + 1]
                    elif tuple(current_loc) in self.bottom_input_state_locations:
                        if self.crossing_signs[self.bottom_input_state_locations.index(tuple(current_loc))] == 1:
                            current_loc = [current_loc[0] - 1, current_loc[1]]
                            up = False
                        else:
                            current_loc = [current_loc[0], current_loc[1] + 1]

                    elif current_loc[1] == self.n_crossings:
                        current_loc = [current_loc[0], 0]
                    else:
                        current_loc = [current_loc[0], current_loc[1] + 1]
                else:
                    if tuple(current_loc) in self.top_output_state_locations:
                        if self.crossing_signs[self.top_output_state_locations.index(tuple(current_loc))] == 1:
                            current_loc = [current_loc[0] + 1, current_loc[1]]
                            up = True
                        else:
                            current_loc = [current_loc[0], current_loc[1] - 1]
                    elif tuple(current_loc) in self.bottom_output_state_locations:
                        if self.crossing_signs[self.bottom_output_state_locations.index(tuple(current_loc))] == 1:
                            current_loc = [current_loc[0] - 1, current_loc[1]]
                            up = True
                        else:
                            current_loc = [current_loc[0], current_loc[1] - 1]
                    elif current_loc[1] == 0:
                        current_loc = [current_loc[0], self.n_crossings]
                    else:
                        current_loc = [current_loc[0], current_loc[1] - 1]
                current_list.append(tuple(old_location))
                unvisited_locations.remove(tuple(old_location))
            component_locations.append(current_list)

        return len(component_locations), component_locations

    def DL(self, L: int) -> bool:
        """Check the DL bound."""
        lower = self.b_resolution_components()[0] - sum([x < 0 for x in self.braid]) - 1
        upper = 1 + sum([x > 0 for x in self.braid]) - self.a_resolution_components()[0]
        return L <= upper and L >= lower

    def load_sign_data(self, sign_data: List[int]):
        """Load sign data into the braid states."""
        from ..constraints.reduction import full_reduce

        for component in range(self.n_components):
            self.strand_signs[component] = sign_data[:self.n_s[component]]
            sign_data = sign_data[self.n_s[component]:]
        self.compute_matrices()
        if self.validate():
            self.generate_position_assignments()
            all_relations = self.get_state_relations()
            relations = full_reduce(all_relations)
            return relations
        return False
