from typing import Optional
from fkcompute.braidstates_links import BraidStates
from fkcompute.relations_links import full_reduce, ilp, print_symbolic_relations
from fkcompute.braids import is_homogeneous_braid

# Import typer
try:
    import typer
except ImportError:
    # Fallback to argparse if typer is not available
    print("Error: typer is required. Please install it with: pip install typer")
    exit(1)

BRAID_TYPE_HOMEGENOUS = 'homogeneous'
BRAID_TYPE_FIBERED = 'fibered'

def braid_type(braid):
    if is_homogeneous_braid(braid):
        return BRAID_TYPE_HOMEGENOUS
    else:
        return BRAID_TYPE_FIBERED

def FK(braid, degree, outfile, inversion_data=None, variables=False):
    braid_type_ = braid_type(braid)
    if braid_type_ == BRAID_TYPE_HOMEGENOUS:
        braid_states = BraidStates(braid)
        all_relations = braid_states.get_state_relations()
        relations = full_reduce(all_relations)
        if variables:
            print_symbolic_relations(degree, relations, braid_states, outfile)
        else:
            solution = ilp(degree, relations, braid_states, outfile)
    elif braid_type_ == BRAID_TYPE_FIBERED:
        if inversion_data is None:
            raise Exception("You must supply inversion data for fibered knots!")
        braid_states = BraidStates(braid)
        braid_states.strand_signs = inversion_data
        braid_states.compute_matrices()
        if braid_states.validate():
            braid_states.generate_position_assignments()
            all_relations = braid_states.get_state_relations()
            relations = full_reduce(all_relations)
            if variables:
                print_symbolic_relations(degree, relations, braid_states, outfile)
            else:
                ilp(degree, relations, braid_states, outfile)
        else:
            raise Exception("The inversion data input doesn't seem to be valid. May you make sure you've input it correctly?")

# -------------------------------------------------------------------------
# Typer App
# -------------------------------------------------------------------------
app = typer.Typer(
    help="A simple CLI tool for FK computation.",
    no_args_is_help=True,
    rich_markup_mode="rich",
)

@app.command()
def main(
    verbose: bool = typer.Option(False, "-v", "--verbose", help='Enable verbose output.'),
    outfile: str = typer.Option(..., "-f", "--outfile", help='File title to write to.'),
    infile: Optional[str] = typer.Option(None, "-i", "--infile", help='Sign assignment input file.'),
    degree: int = typer.Option(..., "-d", "--degree", help='Degree for which sign assignment stays valid.'),
    signs: Optional[str] = typer.Option(None, "-s", "--signs", help='Sign assignment.'),
    partial: Optional[str] = typer.Option(None, "-p", "--partial", help='Sign assignment.'),
    braid: Optional[str] = typer.Option(None, "-b", "--braid", help="The braid."),
    variables: bool = typer.Option(False, "--variables", help='Print symbolic relations instead of computing FK homology.')
) -> None:
    """A simple CLI tool for FK computation."""
    if verbose is None:
        verbose = False
    if infile is None and braid is None:
        raise typer.BadParameter('Need inversion data file to read or a braid!')
    if outfile is None:
        outfile = infile
        
    if infile is not None:
        with open(infile) as file:
            lines = file.readlines()
        sign_assignment = [int(x) for x in lines[1].split(',')]
        if signs is not None:
            sign_assignment = [int(x) for x in signs[1:-1].split(',')]
        braid_list = [int(x) for x in lines[2].split(',')]
    else:
        braid_list = [int(x) for x in braid[1:-1].split(',')]
        if braid_list is None:
            raise typer.BadParameter("Inversion data not found!")
        sign_assignment = None
        
    if partial is None:
        FK(braid_list, degree=degree, outfile=outfile, inversion_data=sign_assignment, variables=variables)
    elif partial == '1':
        # FKlabels function is not defined in the original code, so I'll leave this as a placeholder
        raise typer.BadParameter("FKlabels function is not implemented")
    elif partial == '2':
        FK(braid_list, degree=degree, outfile=outfile, inversion_data=sign_assignment, variables=variables)

if __name__ == "__main__":
    app()