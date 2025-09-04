from fkcompute.braidstates_links import BraidStates
from fkcompute.relations_links import full_reduce, ilp
from fkcompute.braids import is_homogeneous_braid
import argparse

BRAID_TYPE_HOMEGENOUS = 'homogeneous'
BRAID_TYPE_FIBERED = 'fibered'

def braid_type(braid):
    if is_homogeneous_braid(braid):
        return BRAID_TYPE_HOMEGENOUS
    else:
        return BRAID_TYPE_FIBERED

def FK(braid, degree, outfile, inversion_data=None):
    braid_type_ = braid_type(braid)
    if braid_type_ == BRAID_TYPE_HOMEGENOUS:
        braid_states = BraidStates(braid)
        all_relations = braid_states.get_state_relations()
        relations = full_reduce(all_relations)
        solution = ilp(degree, relations, braid_states, outfile)
    elif braid_type_ == BRAID_TYPE_FIBERED:
        if inversion_data == None:
            raise Exception("You must supply inversion data for fibered knots!")
        braid_states = BraidStates(braid)
        braid_states.strand_signs = inversion_data
        braid_states.compute_matrices()
        if braid_states.validate():
            braid_states.generate_position_assignments()
            all_relations = braid_states.get_state_relations()
            relations = full_reduce(all_relations)
            ilp(degree, relations, braid_states, outfile)
        else:
            raise Exception("The inversion data input doesn't seem to be valid. May you make sure you've input it correctly?")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='A simple CLI tool.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output.')
    parser.add_argument('-f', '--outfile', help='File title to write to.')
    parser.add_argument('-i', '--infile', help='Sign assignment input file.')
    parser.add_argument('-d', '--degree', help='Degree for which sign assignment stays valid.')
    parser.add_argument('-s', '--signs', help='Sign assignment.')
    parser.add_argument('-p', '--partial', help='Sign assignment.')
    parser.add_argument('-b', '--braid', help="The braid.")
    args = parser.parse_args()

    if args.outfile == None:
        raise Exception('Need a file title to write to!')

    if args.degree == None:
        raise Exception('Need a degree!')
    else:
        args.degree = int(args.degree)
        
    if args.verbose == None:
        args.verbose = False
    if args.infile == None and args.braid == None:
        raise Exception('Need inversion data file to read or a braid!')
    if args.outfile == None:
        args.outfile = args.infile
    if args.infile is not None:
        with open(args.infile) as file:
            lines = file.readlines()
        sign_assignment = [int(x) for x in lines[1].split(',')]
        if args.signs is not None:
            sign_assignment = [int(x) for x in args.signs[1:-1].split(',')]
        braid = [int(x) for x in lines[2].split(',')]
    else:
        braid = [int(x) for x in args.braid[1:-1].split(',')]
        if braid == None:
            raise Exception("Inversion data not found!")
        sign_assignment = None
    if args.partial == None:
        FK(braid, degree=args.degree, file_title=args.outfile, verbose=args.verbose, inversion_data=sign_assignment)
    elif args.partial == '1':
        FKlabels(braid=braid, degree=args.degree, file_title=args.outfile, inversion_data=sign_assignment)
    elif args.partial == '2':
        FK(braid, degree=args.degree, file_title=args.outfile, verbose=args.verbose, inversion_data=sign_assignment, states_precomputed=True)
