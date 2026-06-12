"""
@file bounded_check.py
@brief Integer Linear Programming (ILP) boundedness verification for Gukov-Manolescu invariant computation

This module provides functionality to determine whether an integer linear programming
problem has bounded feasible solutions. It is used in the context of verifying the
boundedness of linear constraints arising from knot invariant calculations.

The main function constructs a tableau from polynomial constraints and uses Gurobi
optimizer to test boundedness by maximizing each variable in the constraint system.
"""

import gurobipy as gp
from gurobipy import GRB

import numpy as np

from symabelgroup import Symbol

# Configure Gurobi environment to suppress output for cleaner execution
env = gp.Env(empty=True)
env.setParam('OutputFlag', 0)  # Disable solver output messages
env.start()

def integral_bounded(multiples, singlesigns): 
    """
    Check if the integer linear programming problem is bounded.
    
    This function determines whether the ILP problem defined by the given constraints
    has bounded feasible solutions. It constructs a constraint tableau from polynomial
    expressions and tests boundedness by attempting to maximize each variable.
    
    Args:
        multiples: List of polynomial expressions representing constraint coefficients.
                  Each multiple has attributes .var (coefficient array) and methods
                  .is_constant() and .constant() for checking constant terms.
        singlesigns: Dictionary mapping Symbol objects to sign values (+1 or -1).
                    Defines the sign constraints for variables in the system.
    
    Returns:
        bool: True if the ILP problem is bounded (all variables can be maximized
              to finite values), False if unbounded or infeasible.
    
    Algorithm:
        1. Construct tableau matrix from polynomial coefficient arrays
        2. Apply sign transformations based on singlesigns dictionary
        3. Create Gurobi ILP model with non-negativity constraints
        4. Test boundedness by maximizing each variable individually
        5. Return False if any variable is unbounded or if no feasible solution exists
    
    Note: Early termination occurs if any constraint has constant term -1,
          indicating immediate infeasibility.
    """
    n_multiples = len(multiples)
    
    # Get sizes of coefficient arrays for each polynomial constraint
    sizes = [multiple.var.size for multiple in multiples]
    for_fill = max(sizes)  # Maximum array size for padding
    
    # Build constraint tableau matrix
    tableau = []
    for index in range(n_multiples):
        # Early termination: constant -1 indicates infeasible constraint
        if multiples[index].is_constant() and multiples[index].constant() == -1:
            return False
        
        # Pad coefficient arrays to uniform size and add to tableau
        tableau.append(np.concat((multiples[index].var, np.zeros(for_fill - sizes[index]))))
    
    tableau = np.array(tableau)
    
    # Apply sign transformations based on variable sign assignments
    for index in range(1, for_fill):
        if Symbol(index) in singlesigns.keys():
            if singlesigns[Symbol(index)] == -1:
                # For negative-signed variables: negate coefficients and adjust constant term
                tableau[:, index] *= -1
                tableau[:, 0] += tableau[:, index]
    
    # Create Gurobi ILP model for boundedness testing
    model = gp.Model(env=env)
    model.setParam(gp.GRB.Param.PoolSearchMode, 1)  # Enable solution pool mode
    
    # Add integer variables for each symbol in singlesigns
    x = model.addVars(singlesigns.keys(), vtype=GRB.INTEGER)
    
    # Add non-negativity constraints: tableau * x >= 0
    for index in range(n_multiples):
        value = tableau[index][0]  # Constant term
        
        # Add variable terms: sum(coefficient * variable)
        for index_ in range(1, for_fill):
            if tableau[index][index_] != 0:
                value += tableau[index][index_] * x[Symbol(index_)]
        
        # Add constraint: linear expression >= 0
        try:
            model.addConstr(0 <= value)
        except:
            # Skip malformed constraints (defensive programming)
            pass
    
    # Test boundedness by maximizing each variable individually
    for key in singlesigns.keys():
        model.setObjective(x[key], sense=GRB.MAXIMIZE)
        model.optimize()
        
        # Check if optimization was successful (status 2 = optimal solution found)
        if model.Status != 2:
            return False  # Unbounded or infeasible
    
    # Final check: ensure at least one feasible solution exists
    if model.SolCount == 0:
        return False
    
    return True  # Problem is bounded