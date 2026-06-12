"""Generate Mathematica files for FK computation."""

import json
import subprocess
import sys
import tempfile
import os
from typing import Optional, List

try:
    import typer
except ImportError:
    print("Error: typer is required. Install with: pip install typer")
    sys.exit(1)

from fkcompute import BraidStates
from fkcompute.domain.braid.word import is_homogeneous_braid
from fkcompute.domain.constraints.reduction import full_reduce
from fkcompute.solver.ilp import _czech_sign_assignment
from fkcompute.infra.config import load_config_file, parse_int_list


app = typer.Typer(
    help="""
Generate Mathematica files for FK computation.

SIMPLE USAGE:
  mathematica-gen simple "[1,-2,3]" 5
  mathematica-gen simple "[1,-2,3]" 5 --inversion-file inv.json --output out.m

CONFIG FILE USAGE:
  mathematica-gen config config.yaml

CONFIG FILE KEYS:
  braid, degree, inversion (optional), output (optional), run (optional)
""",
    no_args_is_help=True,
)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def make_mathematica_file(braid: list, degree: int, inversion_data: Optional[dict] = None, weight: Optional[int] = None) -> str:
    """Build the Mathematica source string for the given braid and degree."""
    braid_states = BraidStates(braid)

    if inversion_data is not None:
        braid_states.strand_signs = inversion_data
        braid_states.compute_matrices()
        if not braid_states.validate():
            raise ValueError("Invalid inversion data for this braid")
        braid_states.generate_position_assignments()
    elif not is_homogeneous_braid(braid):
        raise ValueError("inversion_data is required for non-homogeneous braids")
    # homogeneous braids: BraidStates already initialised signs in __init__

    all_relations = braid_states.get_state_relations()
    relations = full_reduce(all_relations)

    check = _czech_sign_assignment(degree, relations, braid_states, weight=weight)
    if check is None:
        raise ValueError(f"No valid sign assignment exists at degree {degree}")

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
    for entry in check["multi_var_inequalities"]:
        inequalities += str(entry) + " >= 0, "
    for variable, sign in check["single_var_signs"].items():
        if sign > 0:
            inequalities += str(variable) + " >= 0, "
        else:
            inequalities += str(variable) + " < 0, "
    inequalities = inequalities[:-2] + "};\n\n"

    mathematica_file += inequalities + " res = Collect[Simplify[Normal[Series[Total["
    assignment = check["assignment"]
    c_variables = ["x", "y", "w", "z"]
    R_product = ""
    for i in range(len(braid_states.braid)):
        R_matrix_type = int(braid_states.r_matrices[i][1])
        R_product += "R[" + str(assignment[braid_states.top_input_state_locations[i]]) + ", " \
                         + str(assignment[braid_states.bottom_input_state_locations[i]]) + ", " \
                         + str(assignment[braid_states.top_output_state_locations[i]]) + ", " \
                         + str(assignment[braid_states.bottom_output_state_locations[i]]) + ", " \
                         + str(R_matrix_type) + ", " \
                         + str(c_variables[braid_states.top_crossing_components[i]]) + ", " \
                         + str(c_variables[braid_states.bottom_crossing_components[i]]) + "]*"
    R_product = R_product[:-1]
    v0 = str(c_variables[braid_states.closed_strand_components[0]])
    prefactor = f"({v0}^(1/2)-{v0}^(-1/2))*"
    for c in braid_states.closed_strand_components[1:]:
        prefactor += str(c_variables[c]) + "^(1/2)*"
    if braid_states.n_components == 1:
        prefactor += f"q^({-braid_states.writhe}/4)"
    for k, v in assignment.items():
        if k[1] == 0 and k[0] != 0:
            prefactor += f"q^(-(1+2*({v}))/2)*"
    mathematica_file += prefactor + R_product + "/. # & /@ Solve[ineqs, {"

    variables = list(check["single_var_signs"].keys())
    variables_str = ",".join(str(v) for v in variables)
    mathematica_file += variables_str + "},Integers]]/.{"

    c_variable_inv = ", ".join(f"{v}->{v}^(-1)" for v in c_variables[:braid_states.n_components])
    mathematica_file += c_variable_inv + "} // FunctionExpand,"

    for v in c_variables[:braid_states.n_components]:
        mathematica_file += "{" + str(v) + ",0," + str(degree) + "},"
    mathematica_file = mathematica_file[:-1] + "]], Assumptions -> {q>0,"

    ass = ",".join(f"{v}>0" for v in c_variables[:braid_states.n_components])
    mathematica_file += ass + "}]//Expand,{"

    mathematica_file += ",".join(c_variables[:braid_states.n_components])
    mathematica_file += "}]; \n"
    mathematica_file += "Print[res]"
    return mathematica_file


def run_mathematica_code(code: str) -> str:
    """Execute Mathematica code and return its output."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.m', delete=False) as f:
        f.write(code)
        temp_file = f.name

    try:
        result = subprocess.run(
            ['math', '-script', temp_file],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            return result.stdout.strip()
        return f"Error: {result.stderr.strip()}"
    except subprocess.TimeoutExpired:
        return "Error: Execution timed out"
    except FileNotFoundError:
        return "Error: Mathematica not found. Make sure 'math' is in your PATH"
    finally:
        os.unlink(temp_file)


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

def _emit(content: str, output: Optional[str]) -> None:
    if output:
        with open(output, "w") as f:
            f.write(content)
        typer.echo(f"Written to {output}")
    else:
        typer.echo(content)


def _parse_inversion_file(path: str) -> dict:
    with open(path) as f:
        inv = json.load(f)
    raw = inv.get("inversion_data", inv)
    return {int(k): v for k, v in raw.items()}


def _run_single(config: dict) -> None:
    """Process one computation entry (from CLI args or a config dict)."""
    braid = config.get("braid")
    if not braid:
        typer.echo("Error: 'braid' is required", err=True)
        raise typer.Exit(1)
    degree = config.get("degree")
    if degree is None:
        typer.echo("Error: 'degree' is required", err=True)
        raise typer.Exit(1)

    inversion_data = None
    if config.get("inversion"):
        inversion_data = {int(k): v for k, v in config["inversion"].items()}

    weight = config.get("weight")
    name = config.get("name", "")
    prefix = f"[{name}] " if name else ""

    try:
        code = make_mathematica_file(braid, degree, inversion_data, weight=weight)
    except ValueError as e:
        typer.echo(f"{prefix}Error: {e}", err=True)
        return

    if config.get("run"):
        typer.echo(f"{prefix}{run_mathematica_code(code)}")
    else:
        _emit(code, config.get("output"))


# ---------------------------------------------------------------------------
# CLI commands
# ---------------------------------------------------------------------------

@app.command("simple")
def simple_command(
    braid: str = typer.Argument(..., help='Braid word. Examples: "[1,-2,3]", "1,-2,3"'),
    degree: int = typer.Argument(..., help="Computation degree"),
    inversion_file: Optional[str] = typer.Option(
        None, "--inversion-file", "-i",
        help="Path to inversion JSON file produced by fkcompute",
    ),
    output: Optional[str] = typer.Option(
        None, "--output", "-o",
        help="Output .m file path (default: stdout)",
    ),
    run: bool = typer.Option(False, "--run", help="Execute the generated file with Mathematica"),
    weight: Optional[int] = typer.Option(None, "--weight", "-w", help="Weight parameter for stratified calculation"),
) -> None:
    """Generate a Mathematica file for a single braid."""
    braid_list = parse_int_list(braid)
    if not braid_list:
        raise typer.BadParameter("Could not parse braid into a non-empty list of integers")

    inversion_data = None
    if inversion_file:
        inversion_data = _parse_inversion_file(inversion_file)

    try:
        code = make_mathematica_file(braid_list, degree, inversion_data, weight=weight)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    if run:
        typer.echo(run_mathematica_code(code))
    else:
        _emit(code, output)


@app.command("config")
def config_command(
    config_paths: List[str] = typer.Argument(..., help="Path(s) to YAML or JSON configuration file(s)"),
) -> None:
    """Generate Mathematica file(s) from configuration file(s).

    Supports the same config format as 'fk config', plus two extra keys:
    'output' (path for the .m file) and 'run' (bool, execute with Mathematica).
    Batch mode is supported via a top-level 'computations' list.
    """
    for config_path in config_paths:
        try:
            config_data = load_config_file(config_path)
        except FileNotFoundError as e:
            typer.echo(f"Error: {e}", err=True)
            raise typer.Exit(1)

        computations = config_data.get("computations")
        if computations:
            typer.echo(f"Batch: {len(computations)} computation(s) from {config_path}")
            for comp in computations:
                _run_single(comp)
        else:
            _run_single(config_data)


if __name__ == "__main__":
    app()
