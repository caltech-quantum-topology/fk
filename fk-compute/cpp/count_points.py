"""Count integer points satisfying the Unknown_161658 constraints via Gurobi."""

import sys
sys.path.insert(0, '/home/davep96/work/projects/fk_computation/fk/fk-compute/src')

import gurobipy as gp
from gurobipy import GRB

env = gp.Env(empty=True)
env.setParam('OutputFlag', 0)
env.start()

# Constraints: const + sum_j(coeff_j * x_j) >= 0
# Variables: x[0..9] = n', o, p, r', s, t, u, v, e', i  (all >= 0)

criteria = [
    [2.5, 0, -1,  0, -1,  0,  0,  0,  0, -1,  0],
    [2.5, 0,  0,  2,  1, -1, -1, -1, -1,  0, -1],
]
inequalities = [
    [0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0],
    [0,  0,-1, 1, 0,  0, 0, 0, 0, 0, 0],
    [0,  0, 1, 0, 1,  0, 0, 0, 1,-1,-1],
    [0,  1, 0, 0, 0,  0, 0, 1,-1, 0, 0],
    [0,  0, 0, 1, 0,  0, 0, 0, 0, 0,-1],
    [0,  1, 0, 0, 0,  0, 1, 0,-1, 0, 0],
    [0,  1, 1, 0, 0,  0, 0, 0, 0,-1, 0],
    [0,  0, 1,-1, 0,  0, 0, 0, 1, 0, 0],
    [0, -1, 0, 0, 1,  0, 0, 0, 1, 0, 0],
    [0, -1, 0, 0, 0,  0, 0, 0, 0, 1, 0],
    [0,  1, 1, 0, 0,  0, 0, 0, 0,-1,-1],
    [0,  0, 0, 1, 1,  0, 0, 0, 0,-1,-1],
    [0,  0, 0, 1, 1,  0, 0, 0, 0, 0,-1],
    [0,  1, 0, 0, 0,  1, 0, 0,-1, 0, 0],
]
all_constraints = criteria + inequalities

m = gp.Model(env=env)
m.setParam('PoolSearchMode', 2)   # exhaustive enumeration
m.setParam('PoolSolutions', 10000)  # store up to 10000 solutions

x = m.addVars(10, vtype=GRB.INTEGER, lb=0, name='x')

for c in all_constraints:
    expr = c[0] + gp.quicksum(c[1+j] * x[j] for j in range(10))
    m.addConstr(expr >= 0)

m.setObjective(0)  # feasibility
m.optimize()

if m.SolCount == 0:
    print("No feasible solutions found")
else:
    print(f"Total feasible integer points: {m.SolCount}")
    if m.SolCount <= 30:
        print("Points (n', o, p, r', s, t, u, v, e', i):")
        for k in range(m.SolCount):
            m.setParam('SolutionNumber', k)
            vals = [int(round(x[j].Xn)) for j in range(10)]
            print(f"  {vals}")
