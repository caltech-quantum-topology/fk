from ..braidstates_links_old import BraidStates
import copy # is this bad practice?

# TODO : general resolution components

def a_resolution_components(bs : BraidStates):

    unvisited_locations = copy.deepcopy(bs.state_locations)
    component_locations = []

    while unvisited_locations:
        current_loc = unvisited_locations[0]
        current_list = []
        up = True
        while tuple(current_loc) not in current_list:
            old_location = copy.deepcopy(current_loc)
            if up:
                if tuple(current_loc) in bs.top_input_state_locations:
                    if bs.crossing_signs[bs.top_input_state_locations.index(tuple(current_loc))] == 1:
                        current_loc = [current_loc[0], current_loc[1] + 1]
                    else:
                        current_loc = [current_loc[0] + 1, current_loc[1]]
                        up = False
                elif tuple(current_loc) in bs.bottom_input_state_locations:
                    if bs.crossing_signs[bs.bottom_input_state_locations.index(tuple(current_loc))] == 1:
                        current_loc = [current_loc[0], current_loc[1] + 1]
                    else:
                        current_loc = [current_loc[0] - 1, current_loc[1]]
                        up = False
                elif current_loc[1] == bs.n_crossings:
                    current_loc = [current_loc[0], 0]
                else:
                    current_loc = [current_loc[0], current_loc[1] + 1]
            else:
                if tuple(current_loc) in bs.top_output_state_locations:
                    if bs.crossing_signs[bs.top_output_state_locations.index(tuple(current_loc))] == 1:
                        current_loc = [current_loc[0], current_loc[1] - 1]
                    else:
                        current_loc = [current_loc[0] + 1, current_loc[1]]
                        up = True
                elif tuple(current_loc) in bs.bottom_output_state_locations:
                    if bs.crossing_signs[bs.bottom_output_state_locations.index(tuple(current_loc))] == 1:
                        current_loc = [current_loc[0], current_loc[1] - 1]
                    else:
                        current_loc = [current_loc[0] - 1, current_loc[1]]
                        up = True
                elif current_loc[1] == 0:
                    current_loc = [current_loc[0], bs.n_crossings]
                else:
                    current_loc = [current_loc[0], current_loc[1] - 1]
            current_list.append(tuple(old_location))
            unvisited_locations.remove(tuple(old_location))
        component_locations.append(current_list)
    
    return len(component_locations), component_locations

def b_resolution_components(bs : BraidStates):

    unvisited_locations = copy.deepcopy(bs.state_locations)
    component_locations = []

    while unvisited_locations:
        current_loc = unvisited_locations[0]
        current_list = []
        up = True
        while tuple(current_loc) not in current_list:
            old_location = copy.deepcopy(current_loc)
            if up:
                if tuple(current_loc) in bs.top_input_state_locations:
                    if bs.crossing_signs[bs.top_input_state_locations.index(tuple(current_loc))] == 1:
                        current_loc = [current_loc[0] + 1, current_loc[1]]
                        up = False
                    else:
                        current_loc = [current_loc[0], current_loc[1] + 1]
                elif tuple(current_loc) in bs.bottom_input_state_locations:
                    if bs.crossing_signs[bs.bottom_input_state_locations.index(tuple(current_loc))] == 1:
                        current_loc = [current_loc[0] - 1, current_loc[1]]
                        up = False
                    else:
                        current_loc = [current_loc[0], current_loc[1] + 1]
                    
                elif current_loc[1] == bs.n_crossings:
                    current_loc = [current_loc[0], 0]
                else:
                    current_loc = [current_loc[0], current_loc[1] + 1]
            else:
                if tuple(current_loc) in bs.top_output_state_locations:
                    if bs.crossing_signs[bs.top_output_state_locations.index(tuple(current_loc))] == 1:
                        current_loc = [current_loc[0] + 1, current_loc[1]]
                        up = True
                    else:
                        current_loc = [current_loc[0], current_loc[1] - 1]
                elif tuple(current_loc) in bs.bottom_output_state_locations:
                    if bs.crossing_signs[bs.bottom_output_state_locations.index(tuple(current_loc))] == 1:
                        current_loc = [current_loc[0] - 1, current_loc[1]]
                        up = True
                    else:
                        current_loc = [current_loc[0], current_loc[1] - 1]
                elif current_loc[1] == 0:
                    current_loc = [current_loc[0], bs.n_crossings]
                else:
                    current_loc = [current_loc[0], current_loc[1] - 1]
            current_list.append(tuple(old_location))
            unvisited_locations.remove(tuple(old_location))
        component_locations.append(current_list)
    
    return len(component_locations), component_locations