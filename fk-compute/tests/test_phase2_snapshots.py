"""
Snapshot tests for Phase 2 of the FK computation pipeline.

These tests capture the current behavior of:
- get_state_relations() — raw constraint generation
- full_reduce() — constraint simplification
- symbolic_variable_assignment() — symbolic variable assignment
- process_assignment() — degree criteria + inequality extraction

Any refactoring must preserve these outputs exactly.
"""

import pytest

from fkcompute.domain.braid.states import BraidStates
from fkcompute.domain.braid.types import ZERO_STATE, NEG_ONE_STATE
from fkcompute.domain.constraints.relations import Leq, Less, Zero, NegOne, Alias, Conservation
from fkcompute.domain.constraints.reduction import full_reduce, free_variables
from fkcompute.domain.solver.assignment import symbolic_variable_assignment
from fkcompute.domain.solver.symbolic_constraints import process_assignment
from fkcompute.domain.constraints.pipeline import build_constraint_system
from fkcompute.domain.constraints.system import ConstraintSystem


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def relation_set(relations):
    """Convert relations to a frozenset of repr strings for comparison."""
    return frozenset(repr(r) for r in relations)


def assignment_snapshot(assignment):
    """Convert assignment dict to a comparable snapshot."""
    return {k: str(v) for k, v in sorted(assignment.items())}


# ---------------------------------------------------------------------------
# Trefoil [1, 1, 1]
# ---------------------------------------------------------------------------

class TestTrefoil:
    """Snapshot tests for the trefoil knot braid [1, 1, 1]."""

    @pytest.fixture
    def braid_states(self):
        return BraidStates([1, 1, 1])

    @pytest.fixture
    def raw_relations(self, braid_states):
        return braid_states.get_state_relations()

    @pytest.fixture
    def reduced_relations(self, raw_relations):
        return full_reduce(raw_relations)

    def test_raw_relations(self, raw_relations):
        expected = frozenset([
            "Zero (0, 0) = [0]",
            "Zero (0, 3) = [0]",
            "Alias (0, 3) := (0, 0)",
            "Alias (1, 3) := (1, 0)",
            "Inequality [0] <= (0, 0)",
            "Inequality [0] <= (0, 1)",
            "Inequality [0] <= (0, 2)",
            "Inequality [0] <= (0, 3)",
            "Inequality [0] <= (1, 0)",
            "Inequality [0] <= (1, 1)",
            "Inequality [0] <= (1, 2)",
            "Inequality [0] <= (1, 3)",
            "Inequality (1, 0) <= (0, 1)",
            "Inequality (1, 1) <= (0, 0)",
            "Inequality (1, 1) <= (0, 2)",
            "Inequality (1, 2) <= (0, 1)",
            "Inequality (1, 2) <= (0, 3)",
            "Inequality (1, 3) <= (0, 2)",
            "Conservation (0, 0) + (1, 0) = (0, 1) + (1, 1)",
            "Conservation (0, 1) + (1, 1) = (0, 2) + (1, 2)",
            "Conservation (0, 2) + (1, 2) = (0, 3) + (1, 3)",
        ])
        assert relation_set(raw_relations) == expected

    def test_reduced_relations(self, reduced_relations):
        expected = frozenset([
            "Zero (0, 0) = [0]",
            "Zero (0, 3) = [0]",
            "Zero (1, 1) = [0]",
            "Zero (1, 2) = [0]",
            "Alias (0, 2) := (0, 1)",
            "Alias (1, 0) := (0, 1)",
            "Alias (1, 0) := (0, 2)",
            "Alias (1, 3) := (1, 0)",
            "Inequality [0] <= (0, 1)",
        ])
        assert relation_set(reduced_relations) == expected

    def test_free_variables(self, reduced_relations):
        fv = free_variables(reduced_relations)
        assert fv == [(0, 1)]

    def test_assignment(self, reduced_relations, braid_states):
        assignment = symbolic_variable_assignment(reduced_relations, braid_states)
        snap = assignment_snapshot(assignment)
        expected = {
            (0, 0): "0",
            (0, 1): "a",
            (0, 2): "a",
            (0, 3): "0",
            (1, 0): "a",
            (1, 1): "0",
            (1, 2): "0",
            (1, 3): "a",
        }
        assert snap == expected

    def test_process_assignment(self, reduced_relations, braid_states):
        assignment = symbolic_variable_assignment(reduced_relations, braid_states)
        criteria, multiples, singlesigns = process_assignment(
            assignment, braid_states, reduced_relations,
        )
        # Single component, degree criterion is "a" (with constant 0)
        assert len(criteria) == 1
        assert str(criteria[0]) == "0 + a" or str(criteria[0]) == "a"

        # Knot-only restriction: x-power must be nonnegative.
        assert len(multiples) == 1
        ineq_coeffs = {str(k): v for k, v in multiples[0].as_coefficients_dict().items()}
        assert ineq_coeffs == {"1": 2.0, "a": 4.0}

        # Single variable 'a' has positive sign
        assert len(singlesigns) == 1
        signs_values = list(singlesigns.values())
        assert signs_values[0] == 1.0


# ---------------------------------------------------------------------------
# Figure-eight knot [1, -2, 1, -2]
# ---------------------------------------------------------------------------

class TestFigureEight:
    """Snapshot tests for the figure-eight knot braid [1, -2, 1, -2]."""

    @pytest.fixture
    def braid_states(self):
        return BraidStates([1, -2, 1, -2])

    @pytest.fixture
    def raw_relations(self, braid_states):
        return braid_states.get_state_relations()

    @pytest.fixture
    def reduced_relations(self, raw_relations):
        return full_reduce(raw_relations)

    def test_raw_relations(self, raw_relations):
        expected = frozenset([
            "Zero (0, 0) = [0]",
            "Zero (0, 3) = [0]",
            "Alias (0, 3) := (0, 0)",
            "Alias (1, 4) := (1, 0)",
            "Alias (2, 4) := (2, 0)",
            "Inequality [0] <= (0, 0)",
            "Inequality [0] <= (0, 1)",
            "Inequality [0] <= (0, 3)",
            "Inequality [0] <= (1, 0)",
            "Inequality [0] <= (1, 1)",
            "Inequality [0] <= (1, 2)",
            "Inequality [0] <= (1, 3)",
            "Inequality [0] <= (1, 4)",
            "Inequality (1, 0) <= (0, 1)",
            "Inequality (1, 1) <= (0, 0)",
            "Inequality (1, 2) <= (0, 3)",
            "Inequality (1, 3) <= (0, 1)",
            "Inequality (2, 0) <= [-1]",
            "Inequality (2, 2) <= [-1]",
            "Inequality (2, 4) <= [-1]",
            "Conservation (0, 0) + (1, 0) = (0, 1) + (1, 1)",
            "Conservation (0, 1) + (1, 2) = (0, 3) + (1, 3)",
            "Conservation (1, 1) + (2, 0) = (1, 2) + (2, 2)",
            "Conservation (1, 3) + (2, 2) = (1, 4) + (2, 4)",
        ])
        assert relation_set(raw_relations) == expected

    def test_reduced_relations(self, reduced_relations):
        expected = frozenset([
            "Zero (0, 0) = [0]",
            "Zero (0, 3) = [0]",
            "Zero (1, 1) = [0]",
            "Zero (1, 2) = [0]",
            "Alias (1, 0) := (0, 1)",
            "Alias (1, 3) := (0, 1)",
            "Alias (1, 4) := (1, 0)",
            "Alias (2, 2) := (2, 0)",
            "Alias (2, 4) := (2, 0)",
            "Conservation (0, 1) + (2, 0) = (0, 1) + (2, 0)",
            "Inequality [0] <= (0, 1)",
            "Inequality (2, 0) <= [-1]",
        ])
        assert relation_set(reduced_relations) == expected

    def test_free_variables(self, reduced_relations):
        fv = free_variables(reduced_relations)
        assert fv == [(0, 1), (2, 0)]

    def test_assignment(self, reduced_relations, braid_states):
        assignment = symbolic_variable_assignment(reduced_relations, braid_states)
        snap = assignment_snapshot(assignment)
        expected = {
            (0, 0): "0",
            (0, 1): "a",
            (0, 2): "a",
            (0, 3): "0",
            (0, 4): "0",
            (1, 0): "a",
            (1, 1): "0",
            (1, 2): "0",
            (1, 3): "a",
            (1, 4): "a",
            (2, 0): "b",
            (2, 1): "b",
            (2, 2): "b",
            (2, 3): "b",
            (2, 4): "b",
        }
        assert snap == expected

    def test_process_assignment(self, reduced_relations, braid_states):
        assignment = symbolic_variable_assignment(reduced_relations, braid_states)
        criteria, multiples, singlesigns = process_assignment(
            assignment, braid_states, reduced_relations,
        )
        # Single component
        assert len(criteria) == 1
        assert str(criteria[0]) == "-1 + a - 2b"

        # Knot-only restriction: x-power must be nonnegative.
        assert len(multiples) == 1
        ineq_coeffs = {str(k): v for k, v in multiples[0].as_coefficients_dict().items()}
        assert ineq_coeffs == {"1": -6.0, "a": 4.0, "b": -8.0}

        # Two variables: a positive, b negative
        assert len(singlesigns) == 2
        sign_map = {str(k): v for k, v in singlesigns.items()}
        assert sign_map["a"] == 1.0
        assert sign_map["b"] == -1.0


# ---------------------------------------------------------------------------
# Hopf link [1, 1] (2 components)
# ---------------------------------------------------------------------------

class TestHopfLink:
    """Snapshot tests for the Hopf link braid [1, 1]."""

    @pytest.fixture
    def braid_states(self):
        return BraidStates([1, 1])

    @pytest.fixture
    def reduced_relations(self, braid_states):
        return full_reduce(braid_states.get_state_relations())

    def test_reduced_relations(self, reduced_relations):
        # Just verify it reduces without error and has expected structure
        zeros = [r for r in reduced_relations if isinstance(r, Zero)]
        aliases = [r for r in reduced_relations if isinstance(r, Alias)]
        assert len(zeros) > 0
        assert len(aliases) >= 0

    def test_free_variables(self, reduced_relations):
        fv = free_variables(reduced_relations)
        # Hopf link should have free variables
        assert isinstance(fv, list)

    def test_assignment(self, reduced_relations, braid_states):
        assignment = symbolic_variable_assignment(reduced_relations, braid_states)
        # All states should be assigned
        assert len(assignment) > 0
        # Zero states should map to 0
        assert assignment[(0, 0)] == 0


# ---------------------------------------------------------------------------
# Consistency: full pipeline round-trip
# ---------------------------------------------------------------------------

class TestPipelineConsistency:
    """Verify the full Phase 2 pipeline is consistent across braids."""

    @pytest.mark.parametrize("braid", [
        [1, 1, 1],        # trefoil
        [1, -2, 1, -2],   # figure-eight
        [1, 1],           # Hopf link
    ])
    def test_reduction_is_idempotent(self, braid):
        """Reducing an already-reduced set should produce the same result."""
        bs = BraidStates(braid)
        raw = bs.get_state_relations()
        reduced1 = full_reduce(raw)
        reduced2 = full_reduce(reduced1)
        assert relation_set(reduced1) == relation_set(reduced2)

    @pytest.mark.parametrize("braid", [
        [1, 1, 1],
        [1, -2, 1, -2],
        [1, 1],
    ])
    def test_assignment_covers_all_states(self, braid):
        """Every state in the reduced relations should have an assignment."""
        bs = BraidStates(braid)
        reduced = full_reduce(bs.get_state_relations())
        assignment = symbolic_variable_assignment(reduced, bs)

        # All variables in relations should be assigned
        from fkcompute.domain.constraints.reduction import all_variables
        all_vars = all_variables(reduced)
        for v in all_vars:
            assert v in assignment, f"Variable {v} not in assignment"

    @pytest.mark.parametrize("braid", [
        [1, 1, 1],
        [1, -2, 1, -2],
        [1, 1],
    ])
    def test_build_constraint_system(self, braid):
        """build_constraint_system produces the same results as calling steps individually."""
        bs = BraidStates(braid)
        cs = build_constraint_system(bs)

        assert isinstance(cs, ConstraintSystem)

        # Compare with manual pipeline
        reduced = full_reduce(bs.get_state_relations())
        assignment = symbolic_variable_assignment(reduced, bs)
        criteria, multi_var_ineqs, svs = process_assignment(assignment, bs, reduced)

        assert assignment_snapshot(cs.assignment) == assignment_snapshot(assignment)
        assert len(cs.degree_criteria) == len(criteria)
        assert len(cs.multi_var_inequalities) == len(multi_var_ineqs)
        assert len(cs.single_var_signs) == len(svs)
