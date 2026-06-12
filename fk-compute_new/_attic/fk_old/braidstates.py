"""
BraidStates module for knot theory computations in Gukov-Manolescu invariants.

This module defines the BraidStates class which represents the combined data of:
* A braid representation (sequence of generators)
* The states of the braid at each position
* Known relations among states for constraint-based computation

The class handles the geometric and algebraic structure of braids, including:
- Strand tracking and component identification
- Crossing matrix computation for R-matrix types
- State relation generation for algebraic constraints
- Sign assignment validation and position mapping

Used as input to the Gukov-Manolescu invariant computation pipeline.
"""

from relations import Zero, Nunity, Alias, Conservation, Leq, free_variables, full_reduce
import itertools
from stateliteral import NUNITY_STATE, ZERO_STATE

class BraidStates:
    """
    Represents a braid with its states and geometric/algebraic structure.
    
    This class analyzes a braid representation and computes all necessary
    geometric and algebraic data for Gukov-Manolescu invariant computation,
    including strand tracking, component identification, and state relations.
    """
    
    def __init__(self, braid):
        """
        Initialize BraidStates from a braid word.
        
        Args:
            braid: List of integers representing braid generators.
                  Positive integers represent positive crossings,
                  negative integers represent negative crossings.
                  Generator i acts on strands i-1 and i (0-indexed).
        """
        self.braid = braid
        self.max_strand = max(abs(x) for x in braid)  # Highest strand index in braid
        self.strands = list(range(0, self.max_strand + 1))  # All strand indices [0, 1, ..., max]
        self.n_strands = self.max_strand + 1  # Total number of strands
        self.n_crossings = len(braid)  # Number of crossings in braid
        self.braid_group_generators = list(range(1, self.max_strand + 1))  # Valid generator indices
        
        # Extract crossing signs from braid word
        self.crossing_signs = []
        for g in braid:
            if g < 0:
                self.crossing_signs.append(-1)  # Negative crossing
            else:
                self.crossing_signs.append(+1)  # Positive crossing
        self.writhe = sum(self.crossing_signs)  # Total signed crossing count

        # Grid locations for states: (strand_index, crossing_index)
        # States exist at each strand position before and after each crossing
        self.state_locations = [(i, j) for i in self.strands for j in range(0, self.n_crossings + 1)]
        
        # Crossing input/output locations for geometric analysis
        # Top input: upper strand entering each crossing (generator - 1 for 0-indexing)
        self.top_input_state_locations = [(abs(braid[index]) - 1, index) for index in range(self.n_crossings)]
        # Bottom input: lower strand entering each crossing (top + 1)
        self.bottom_input_state_locations = [(x + 1, y) for (x, y) in self.top_input_state_locations]
        # Output locations: same strands after crossing (position + 1)
        self.top_output_state_locations = [(x, y + 1) for (x, y) in self.top_input_state_locations]
        self.bottom_output_state_locations = [(x, y + 1) for (x, y) in self.bottom_input_state_locations]
        
        # Matrix indexing for crossing data: [output][input] positions
        # loc[crossing][output_strand][input_strand] = (strand, position)
        self.loc = [[[(x, y + 1), (x + 1, y + 1)], [(x, y), (x + 1, y)]] for (x, y) in self.top_input_state_locations]

        # State tracking data structures
        state_info = {}  # Metadata for each state
        state_at_location = {}  # Map from (strand, position) to state identifier
        
        # Initialize states for each strand
        for strand in self.strands:
            state = (strand, 0)  # Initial state at position 0
            state_info[state] = {}
            state_at_location[(strand, 0)] = (strand, 0)
            
            # Create states only at positions where strand participates in crossings
            for j in range(1, self.n_crossings + 1):
                # Check if current crossing involves this strand
                if abs(braid[j-1]) == strand or abs(braid[j-1]) == strand + 1:
                    state = (strand, j)  # New state after crossing
                    state_info[state] = {}
                # Map location to current state (unchanged if no crossing involvement)
                state_at_location[(strand, j)] = state
        
        self.state_info = state_info
        self._state_at_location = state_at_location
        self.states = list(sorted(set(state_at_location.values())))  # Unique states

        # Build equivalence classes: group locations that map to same state
        self.state_equivalence_classes = {}
        for state in self.states:
            self.state_equivalence_classes[state] = []
        for key in state_at_location.keys():
            self.state_equivalence_classes[state_at_location[key]].append(key)

        # Component analysis: trace connected strand paths through the braid
        # Each component is a connected path that may close into a loop
        self.component_locations = []  # List of location lists, one per component
        self.strand_locations = dict()  # Map component -> list of segments with endpoints
        index = 0
        component = 0
        
        # Trace each component by following strand connections through crossings
        while index < self.n_strands:
            current_loc = [index, 0]  # Start at strand 'index', position 0
            # Skip if this location is already part of a traced component
            if tuple(current_loc) not in list(itertools.chain.from_iterable(self.component_locations)):
                done = False
                prefix = []  # Segment endpoint type from previous segment
                self.strand_locations[component] = []  # Segments for this component
                current_list = []  # All locations in this component
                current_list_strand = []  # Current segment being traced
                
                # Trace the component path until we return to a visited location
                while not done:
                    current_list.append(tuple(current_loc))
                    current_list_strand.append(tuple(current_loc))
                    
                    # Check crossing participation and determine next location
                    if tuple(current_loc) in self.bottom_input_state_locations:
                        # Bottom strand at crossing: goes to top strand of next position
                        self.strand_locations[component].append(prefix + current_list_strand + ['R'])
                        current_loc = [current_loc[0] - 1, current_loc[1] + 1]
                        current_list_strand = []
                        prefix = ['L']
                    elif tuple(current_loc) in self.top_input_state_locations:
                        # Top strand at crossing: goes to bottom strand of next position
                        self.strand_locations[component].append(prefix + current_list_strand + ['L'])
                        current_loc = [current_loc[0] + 1, current_loc[1] + 1]
                        current_list_strand = []
                        prefix = ['R']
                    elif current_loc[1] == self.n_crossings:
                        # Reached end: wrap back to beginning of same strand
                        current_loc = [current_loc[0], 0]
                    else:
                        # No crossing: continue along same strand
                        current_loc = [current_loc[0], current_loc[1] + 1]
                    
                    # Check for cycle completion
                    if tuple(current_loc) in current_list:
                        # Close the loop by connecting final segment to first
                        self.strand_locations[component][0] = prefix + current_list_strand + self.strand_locations[component][0]
                        done = True
                
                # Store component if non-empty
                if current_list != []:
                    self.component_locations.append(current_list)
                    component += 1
            index += 1

        # Component mapping and analysis
        self.component_locations_dict = {index: list_ for (index, list_) in enumerate(self.component_locations)}
        self.n_components = len(self.component_locations)  # Total number of connected components
        
        # Component identification for crossing analysis
        self.top_crossing_components = [self.get_component(location) for location in self.top_input_state_locations]
        self.bottom_crossing_components = [self.get_component(location) for location in self.bottom_input_state_locations]
        self.closed_strand_components = [self.get_component(location) for location in [(index, 0) for index in range(self.n_strands)]]

        # Extract and process strand segment data
        # strand_types: endpoint types (L/R) indicating crossing participation
        self.strand_types = {key: [[x[0], x[-1]] for x in value] for (key, value) in self.strand_locations.items()}
        # Remove endpoint markers, keep only location sequences
        self.strand_locations = {key: [x[1:-1] for x in value] for (key, value) in self.strand_locations.items()}
        # Extract start/end locations for each segment
        self.strand_endpoints = {key: [[x[0], x[-1]] for x in value] for (key, value) in self.strand_locations.items()}
        # Convert locations to crossing indices for matrix computation
        self.endpoint_crossing_indices = {key: [[x[0][1] - 1, x[-1][1]] for x in value] for (key, value) in self.strand_endpoints.items()}
        
        # Segment counting
        self.n_s = {key: len(value) for (key, value) in self.strand_locations.items()}  # Segments per component
        self.n_s_total = sum(list(self.n_s.values()))  # Total segments = 2 * n_crossings + n_components
        
        # Initialize strand signs (will be set by load_sign_data)
        self.strand_signs = {c: [1 for _ in range(self.n_s[c])] for c in range(self.n_components)}
        
    def compute_r_matrices(self):
        """
        Determine R-matrix type for each crossing based on crossing signs and matrices.
        
        R-matrix classification follows Kim's prescription:
        - Positive crossings (sign +1): R1 (both signs same), R2 (mixed signs)
        - Negative crossings (sign -1): R3 (mixed signs), R4 (both signs same)
        
        Returns:
            0 if all crossings have valid R-matrix types, None if invalid
        """
        self.r_matrices = []
        for index in range(self.n_crossings):
            crossing_sign = self.crossing_signs[index]
            
            if crossing_sign == 1:  # Positive crossing
                if self.matrices[index][0][1] * self.matrices[index][1][0] > 0:
                    self.r_matrices.append('R1')  # Same signs
                elif self.matrices[index][0][1] == 1 and self.matrices[index][1][0] == -1:
                    self.r_matrices.append('R2')  # Mixed signs (1, -1)
                else:
                    return None  # Invalid configuration
                    
            elif crossing_sign == -1:  # Negative crossing
                if self.matrices[index][0][0] * self.matrices[index][1][1] > 0:
                    self.r_matrices.append('R4')  # Same signs
                elif self.matrices[index][0][0] == 1 and self.matrices[index][1][1] == -1:
                    self.r_matrices.append('R3')  # Mixed signs (1, -1)
                else:
                    return None  # Invalid configuration
        return 0  
    
    def validate(self):
        """
        Validate crossing matrices for consistency.
        
        Checks two conditions:
        1. Each crossing preserves strand count (input 1s = output 1s)
        2. All crossings have valid R-matrix classifications
        
        Returns:
            True if valid, False otherwise
        """
        for index in range(self.n_crossings):
            # Check strand conservation: input count = output count
            if self.matrices[index][0].count(1) != self.matrices[index][1].count(1):
                return False
        
        # Check R-matrix validity
        if self.compute_r_matrices() is not None:
            return True
        return False

    def get_component(self, location):
        """
        Find which component contains the given location.
        
        Args:
            location: Tuple (strand, position) to locate
            
        Returns:
            Component index containing the location, or None if not found
        """
        for (key, value) in self.component_locations_dict.items():
            if location in value:
                return key

    def sign_at_position(self, position):
        """
        Get the sign assignment at a specific position from crossing matrices.
        
        Args:
            position: Tuple (strand, crossing) position
            
        Returns:
            Sign value (1 or -1) at the position, or None if not found
        """
        for index in range(self.n_crossings):
            for index1 in range(2):  # Output strands
                for index2 in range(2):  # Input strands
                    if self.loc[index][index1][index2] == position:
                        return self.matrices[index][index1][index2]

    def generate_position_assignments(self):
        """
        Generate sign assignments for all positions based on segment signs.
        
        Each position inherits the sign from the segment that starts at that position.
        This creates a complete mapping from (strand, position) to sign values.
        """
        self.sign_assignment = dict()
        for component in range(self.n_components):
            for index in range(self.n_s[component]):  # For each segment in component
                # Get sign from first position of segment
                sign = self.sign_at_position(self.strand_locations[component][index][0])
                # Assign this sign to all positions in the segment
                for index_ in range(len(self.strand_locations[component][index])):
                    self.sign_assignment[self.strand_locations[component][index][index_]] = sign

    def compute_matrices(self):
        """
        Compute crossing matrices from strand signs and endpoint types.
        
        Each crossing has a 2x2 matrix representing strand flow:
        - matrices[i][0] = output strand values (top/bottom after crossing)
        - matrices[i][1] = input strand values (top/bottom before crossing)
        
        The matrix entries are populated from segment signs at crossing endpoints.
        """
        # Initialize matrices: [crossing][output/input][top/bottom]
        self.matrices = [[[None, None], [None, None]] for _ in range(self.n_crossings)]
        
        for component in range(self.n_components):
            for index in range(self.n_s[component]):  # For each segment
                # Determine input side based on segment start type
                if self.strand_types[component][index][0] == 'L':
                    bottom = 0  # Left side corresponds to bottom input
                elif self.strand_types[component][index][0] == 'R':
                    bottom = 1  # Right side corresponds to top input
                    
                # Determine output side based on segment end type
                if self.strand_types[component][index][1] == 'L':
                    top = 0  # Left side corresponds to bottom output
                elif self.strand_types[component][index][1] == 'R':
                    top = 1  # Right side corresponds to top output
                
                # Set matrix entries for start and end crossings of this segment
                start_crossing = self.endpoint_crossing_indices[component][index][0]
                end_crossing = self.endpoint_crossing_indices[component][index][1]
                sign = self.strand_signs[component][index]
                
                self.matrices[start_crossing][0][bottom] = sign  # Output at start
                self.matrices[end_crossing][1][top] = sign       # Input at end

    def get_state(self, location):
        """
        Get the state identifier at a given location.
        
        Args:
            location: Tuple (strand, position)
            
        Returns:
            State identifier at the location, or None if not found
        """
        return self._state_at_location.get(location)

    def get_state_relations(self):
        """
        Generate all algebraic relations between states for constraint solving.
        
        Creates relations of different types:
        1. Boundary conditions: Zero/Nunity states at strand endpoints
        2. Closure conditions: Alias relations connecting braid ends
        3. Sign constraints: Leq relations based on position sign assignments
        4. Crossing relations: R-matrix constraints and conservation laws
        
        Returns:
            List of relation objects for use in algebraic constraint solving
        """
        relations = []
        
        # Boundary conditions based on first segment sign
        if self.strand_signs[0][0] == 1:
            relations.append(Zero(self.get_state((0, 0))))
            relations.append(Zero(self.get_state((0, len(self.braid)))))
        else:
            relations.append(Nunity(self.get_state((0, 0))))
            relations.append(Nunity(self.get_state((0, len(self.braid)))))
    
        # Closure conditions: connect braid endpoints (periodic boundary)
        for i in self.strands:
            relations.append(Alias(self.get_state((i, len(self.braid))), self.get_state((i, 0))))
                             
        # Sign-based constraints for each state position
        for (i, j) in self.state_info.keys():
            if self.sign_assignment[(i, j)] == 1:
                relations.append(Leq(ZERO_STATE, (i, j)))  # 0 ≤ state
            elif self.sign_assignment[(i, j)] == -1:
                relations.append(Leq((i, j), NUNITY_STATE))  # state ≤ -1
            else:
                raise ValueError(f"The sign at position {(i, j)} is not 1 or -1! This isn't supposed to be the case!")
        
        # Crossing constraints based on R-matrix types
        for j, gen in enumerate(self.braid):
            gen_abs = abs(gen)
            
            if self.r_matrices[j] == 'R1':
                # R1 type relations (positive crossing, same signs)
                relations.append(Leq(self.get_state((gen_abs, j+1)), self.get_state((gen_abs-1, j))))
                relations.append(Leq(self.get_state((gen_abs, j)), self.get_state((gen_abs-1, j+1))))
            elif self.r_matrices[j] == 'R4':
                # R4 type relations (negative crossing, same signs)
                relations.append(Leq(self.get_state((gen_abs-1, j)), self.get_state((gen_abs, j+1))))
                relations.append(Leq(self.get_state((gen_abs-1, j+1)), self.get_state((gen_abs, j))))
    
            # Conservation relation: input state sum = output state sum
            relations.append(Conservation(
                [self.get_state((gen_abs-1, j)), self.get_state((gen_abs, j))],      # Input states
                [self.get_state((gen_abs-1, j+1)), self.get_state((gen_abs, j+1))]  # Output states
            ))     
        return relations

    def reduced_relations(self):
        """
        Get algebraically reduced state relations.
        
        Returns:
            List of reduced relation objects after elimination of redundancies
        """
        return full_reduce(self.get_state_relations())

    def free_variables(self):
        """
        Find free variables in the reduced relation system.
        
        Returns:
            List of state variables that remain free after constraint reduction,
            excluding special boundary states (ZERO_STATE, NUNITY_STATE)
        """
        return [x for x in free_variables(self.reduced_relations()) if x != ZERO_STATE and x != NUNITY_STATE]
    
    def load_sign_data(self, sign_data):
        """
        Load strand sign assignments and validate the resulting configuration.
        
        Args:
            sign_data: Flat list of sign values (+1/-1) for all segments
                      ordered by component, then by segment within component
        
        Returns:
            Reduced relations list if valid configuration, False if invalid
            
        Process:
            1. Distribute sign data to components and segments
            2. Compute crossing matrices from signs and geometry
            3. Validate matrix consistency and R-matrix types
            4. Generate position assignments and state relations
            5. Return reduced relation system for constraint solving
        """
        # Distribute sign data across components
        for component in range(self.n_components):
            self.strand_signs[component] = sign_data[:self.n_s[component]]
            sign_data = sign_data[self.n_s[component]:]
        
        # Compute and validate matrices
        self.compute_matrices()
        if self.validate():
            self.generate_position_assignments()
            all_relations = self.get_state_relations()
            relations = full_reduce(all_relations)
            return relations
        return False
    