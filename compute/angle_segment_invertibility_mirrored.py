from relations import symbolic_variable_assignment, process_assignment, angle_variables, angle_segment_invertibility
from relations import full_reduce
from braidstates import BraidStates
from braids import is_homogeneous_braid

def angle_segment_invertibility_wrapper(relations, braid_states):

    assignment = symbolic_variable_assignment(relations, braid_states)
    _, _, singlesigns = process_assignment(assignment, braid_states, relations)
    assignm, mapping = angle_variables(braid_states, assignment)
    return angle_segment_invertibility(
        singlesigns,
        assignm,
        mapping,
        braid_states,
        assignment
    )

BRAID_TYPE_HOMEGENOUS = 'homogeneous'
BRAID_TYPE_FIBERED = 'fibered'

def braid_type(braid):
    if is_homogeneous_braid(braid):
        return BRAID_TYPE_HOMEGENOUS
    else:
        return BRAID_TYPE_FIBERED

def call_wrapper(braid, inversion_data=None):

    braid_type_ = braid_type(braid)
    if braid_type_ == BRAID_TYPE_HOMEGENOUS:
        braid_states = BraidStates(braid)
        all_relations = braid_states.get_state_relations()
        relations = full_reduce(all_relations)
        a, b = angle_segment_invertibility_wrapper(relations, braid_states)
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
            a, b = angle_segment_invertibility_wrapper(relations, braid_states)
        else:
            raise Exception("The inversion data input doesn't seem to be valid. May you make sure you've input it correctly?")

    return a, b

def writer(element, outfile):
    a, b = call_wrapper(
        element["braid"], 
        element["inversion_data"]
    )
    with open(outfile + ".csv", 'w') as file:
        for row in a:
            for e in row:
                file.write(str(int(e)) + ",")
            file.write("\n")
    with open(outfile + "_inv.csv", 'w') as file:
        for row in b:
            for e in row:
                file.write(str(int(e)) + ",")
            file.write("\n")

if __name__ == "__main__":

    import json
    infile = "data/fk_fibered_inversion_data_mirrored.json"

    data = json.load(
        open(infile, "r")
    )

    failure_count = 0
    for element in list(data.values()):
        try:
            writer(
                element, 
                "angle_segment_invertibility_mirrored/" + element["knotinfo_id"]
            )
        except:
            failure_count += 1
    print("failure count:", failure_count)