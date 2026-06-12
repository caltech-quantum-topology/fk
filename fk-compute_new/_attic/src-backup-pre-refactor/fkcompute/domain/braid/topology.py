"""
BraidTopology class for pure topological braid data.

This class represents the topological structure of a braid closure,
independent of any sign assignment. It contains:
- Braid word and generator information
- State locations and equivalence classes
- Component structure (which strands form which components)
- Strand endpoint and type information
"""

import string
import copy
import itertools
from typing import List, Dict, Optional, Tuple, Any


class BraidTopology:
    """
    Pure topological data for a braid closure.

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

    def get_component(self, location: Tuple[int, int]) -> Optional[int]:
        """Get the component index for a given location."""
        for (key, value) in self.component_locations_dict.items():
            if location in value:
                return key
        return None

    def get_state(self, location: Tuple[int, int]) -> Optional[Tuple[int, int]]:
        """Get the state at a given location."""
        return self._state_at_location.get(location)

    def label_from_location(self, location) -> str:
        """Get the label for a location."""
        if type(location) != tuple:
            return str(location)
        return self.state_info[self.get_state(location)].get('label') or str(location)

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

    def Lobb(self, L: int) -> bool:
        """Check the Lobb bound."""
        lower = sum([-x not in self.braid for x in range(1, self.n_strands)])
        upper = sum([x not in self.braid for x in range(1, self.n_strands)])
        U = self.writhe + self.n_strands + 1 - 2 * lower
        Delta = self.n_strands + 1 - upper - lower
        return L <= U and L >= U - 2 * Delta

    def DL(self, L: int) -> bool:
        """Check the DL bound."""
        lower = self.b_resolution_components()[0] - sum([x < 0 for x in self.braid]) - 1
        upper = 1 + sum([x > 0 for x in self.braid]) - self.a_resolution_components()[0]
        return L <= upper and L >= lower
