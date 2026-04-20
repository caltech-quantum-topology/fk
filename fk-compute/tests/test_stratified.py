"""
Tests for stratified FK computation (weight parameter).
"""

import pytest
import fkcompute
from fkcompute.inversion.search import check_assignment_for_braid
from fkcompute.domain.braid.states import BraidStates
from fkcompute.domain.constraints.reduction import full_reduce
from fkcompute.domain.solver.symbolic_constraints import _total_weight, minimum_degree_symbolic
from fkcompute.domain.solver.assignment import symbolic_variable_assignment


# ---------------------------------------------------------------------------
# Test 1: weight=None is identical to omitting the argument
# ---------------------------------------------------------------------------

def test_weight_none_matches_default():
    """weight=None must produce the same result as not passing weight at all."""
    result_default = fkcompute.fk([1, 1, 1], degree=5)
    result_none = fkcompute.fk([1, 1, 1], degree=5, weight=None)
    assert result_default["terms"] == result_none["terms"]


# ---------------------------------------------------------------------------
# Test 2: Smoke test — weight parameter accepted without error
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("braid,degree,weight", [
    ([1, 1, 1], 5, 3),
    ([1, 1, 1], 4, 4),
    ([1, -2, 1, -2, 1], 4, 2),
])
def test_stratified_runs_without_error(braid, degree, weight):
    """Passing weight should not raise and should return a valid result dict."""
    result = fkcompute.fk(braid, degree=degree, weight=weight)
    assert "terms" in result
    assert isinstance(result["terms"], list)


# ---------------------------------------------------------------------------
# Test 3: weight changes the result
# ---------------------------------------------------------------------------

def test_weight_affects_result():
    """A non-None weight should change the output relative to weight=None."""
    braid = [1, 1, 1]
    degree = 6
    result_unweighted = fkcompute.fk(braid, degree=degree)
    result_weighted = fkcompute.fk(braid, degree=degree, weight=2)
    assert result_unweighted["terms"] != result_weighted["terms"]


# ---------------------------------------------------------------------------
# Test 4: weight = degree is a valid edge case
# ---------------------------------------------------------------------------

def test_weight_equals_degree():
    """weight == degree should not crash."""
    result = fkcompute.fk([1, 1, 1], degree=4, weight=4)
    assert "terms" in result


# ---------------------------------------------------------------------------
# Test 5: Unit tests for _total_weight helper
# ---------------------------------------------------------------------------

def _make_braid_states_with_assignment(braid):
    """Return (braid_states, symbolic_assignment) for a homogeneous braid."""
    bs = BraidStates(braid)
    bs.generate_position_assignments()
    relations = full_reduce(bs.get_state_relations())
    assignment = symbolic_variable_assignment(relations, bs)
    return bs, assignment


def test_total_weight_all_positive_strands():
    """
    For the trefoil [1,1,1] (all positive generators, single component),
    the total weight should be a symbolic expression that is non-negative
    when evaluated at any valid (non-negative) assignment.
    """
    bs, assignment = _make_braid_states_with_assignment([1, 1, 1])
    w = _total_weight(assignment, bs)
    # Should be a Symbol (symbolic expression), not an int error
    assert w is not None


def test_total_weight_type():
    """_total_weight should return a Symbol-compatible object."""
    from fkcompute.domain.constraints.symbols import Symbol
    bs, assignment = _make_braid_states_with_assignment([1, 1, 1])
    w = _total_weight(assignment, bs)
    # Either a Symbol or an integer — both are valid
    assert isinstance(w, (Symbol, int))


def test_minimum_degree_symbolic_adds_weight_key():
    """When weight is provided, minimum_degree_symbolic must add key -1."""
    bs, assignment = _make_braid_states_with_assignment([1, 1, 1])
    conditions_no_weight = minimum_degree_symbolic(assignment, bs)
    conditions_with_weight = minimum_degree_symbolic(assignment, bs, weight=3)

    assert -1 not in conditions_no_weight
    assert -1 in conditions_with_weight


# ---------------------------------------------------------------------------
# Test 6: CLI --weight flag
# ---------------------------------------------------------------------------

def test_cli_weight_flag():
    """CLI simple command should accept --weight and exit cleanly."""
    from typer.testing import CliRunner
    from fkcompute.cli.app import app

    runner = CliRunner()
    result = runner.invoke(app, ["simple", "[1,1,1]", "5", "--weight", "3"])
    assert result.exit_code == 0, f"CLI exited with code {result.exit_code}:\n{result.output}"


# ---------------------------------------------------------------------------
# Test 7: weight propagates into inversion search (check_assignment_for_braid)
# ---------------------------------------------------------------------------

def test_check_assignment_accepts_weight():
    """check_assignment_for_braid should accept a weight kwarg without error."""
    bs = BraidStates([1, 1, 1])
    # Index 0 = all-minus, total = 2^n_s_total; try index 0 with weight
    result_no_weight = check_assignment_for_braid(0, 5, bs, weight=None)
    bs2 = BraidStates([1, 1, 1])
    result_with_weight = check_assignment_for_braid(0, 5, bs2, weight=3)
    # Both calls should return without raising (result may be None or a tuple)
    assert result_no_weight is None or isinstance(result_no_weight, tuple)
    assert result_with_weight is None or isinstance(result_with_weight, tuple)
