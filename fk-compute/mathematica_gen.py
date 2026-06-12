import pandas as pd
from fkcompute.braidstates_links import BraidStates
from fkcompute.relations_links import full_reduce,czech_sign_assignment
from fkcompute.relations_links import Symbol
import subprocess
import tempfile
import os
import fkcompute
import yaml
from tqdm import tqdm

def make_mathematica_file(braid, degree, inversion_data):
    braid_states = BraidStates(braid)
    braid_states.strand_signs = inversion_data
    braid_states.compute_matrices()
    if braid_states.validate():
        braid_states.generate_position_assignments()
        all_relations = braid_states.get_state_relations()
        relations = full_reduce(all_relations)
    else:
        raise ValueError("Invalid inversion data computed for this braid")
    check = czech_sign_assignment(degree,relations,braid_states)
    mathematica_file = '''
    << QAlg``
    (*IntegerSolutions:list integer points satisfying inequalities.Usage:\
    IntegerSolutions[ineqs,{x,y,...}] IntegerSolutions[ineqs,{x,y,...},\
    "Bounds"->{{xmin,xmax},{ymin,ymax},...}] IntegerSolutions[ineqs,{x,y,...},\
    "Bounds"->{{x,xmin,xmax},{y,ymin,ymax},...}]*)
    ClearAll[IntegerSolutions];
    Options[IntegerSolutions] = {"Bounds" -> Automatic};
    '''
    
    inequalities = "ineqs = {"
    for value in check["criteria"].values():
        inequalities += str(value) + " >= 0, "
    for entry in check["multiples"]:
        inequalities += str(entry) + " >= 0, "
    for variable, sign in check["single_signs"].items():
        if sign>0:
            inequalities += str(variable) + " >= 0, " 
        else:
            inequalities += str(variable) + " < 0, "
    inequalities = inequalities[:-2] + "};\n\n"
    
    mathematica_file += inequalities + " res = Collect[Simplify[Normal[Series[Total["
    assignment = check["assignment"]
    c_variables = ["x","y","w","z"]
    R_product = ""
    seen = set() 
    for i in range(len(braid_states.braid)):
        R_matrix_type = int(braid_states.r_matrices[i][1])
        R_product += "R["+str(assignment[braid_states.top_input_state_locations[i]])+", "\
                         +str(assignment[braid_states.bottom_input_state_locations[i]])+", "\
                         +str(assignment[braid_states.top_output_state_locations[i]])+", "\
                         +str( assignment[braid_states.bottom_output_state_locations[i]])+", "\
                         +str(R_matrix_type)+", "\
                         +str(c_variables[braid_states.top_crossing_components[i]])+", "\
                         +str(c_variables[braid_states.bottom_crossing_components[i]])+"]*"
    R_product = R_product[:-1]
    v0 = str(c_variables[braid_states.closed_strand_components[0]])
    prefactor = f"({v0}^(1/2)-{v0}^(-1/2))*"
    for c in braid_states.closed_strand_components[1:]:
        prefactor += str(c_variables[c]) + "^(1/2)*"
    if braid_states.n_components == 1:
        prefactor += f"q^({-braid_states.writhe}/4)"
    for k,v in assignment.items():
        if k[1] == 0 and k[0] != 0:
            prefactor += f"q^(-(1+2*({v}))/2)*"
    mathematica_file += prefactor + R_product + "/. # & /@ Solve[ineqs, {"
    
    variables = list(check["single_signs"].keys())
    variables_str = ""
    for v in variables:
        variables_str += str(v)+","
    variables_str = variables_str[:-1]
    mathematica_file += variables_str + "},Integers]]/.{"
    
    c_variable_inv = ""
    for v in c_variables[:braid_states.n_components]:
        c_variable_inv += str(v) + "->" + str(v) + "^(-1), "
    c_variable_inv = c_variable_inv[:-2] + "}"
    
    mathematica_file += c_variable_inv + "// FunctionExpand,"
    for v in c_variables[:braid_states.n_components]:
        v_expr_str = "{" + str(v) + ",0," + str(degree) + "},"
        mathematica_file += v_expr_str
    mathematica_file = mathematica_file[:-1] + "]], Assumptions -> {q>0,"
    
    ass = ""
    for v in c_variables[:braid_states.n_components]:
        ass += str(v) + "> 0,"
    ass = ass[:-1] + "}]//Expand,{"
    mathematica_file += ass
    
    for v in c_variables[:braid_states.n_components]:
        mathematica_file += str(v) + ","
    mathematica_file = mathematica_file[:-1] + "}]; \n"
    mathematica_file += "Print[res]"
    return mathematica_file

def run_mathematica_code(code):
    """
    Execute Mathematica code and return the result.

    Args:
        code (str): Mathematica code to execute

    Returns:
        str: Output from Mathematica execution
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.m', delete=False) as f:
        f.write(code)
        temp_file = f.name

    try:
        result = subprocess.run(
            ['math', '-script', temp_file],
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            return result.stdout.strip()
        else:
            return f"Error: {result.stderr.strip()}"

    except subprocess.TimeoutExpired:
        return "Error: Execution timed out"
    except FileNotFoundError:
        return "Error: Mathematica not found. Make sure 'math' is in your PATH"
    finally:
        os.unlink(temp_file)
