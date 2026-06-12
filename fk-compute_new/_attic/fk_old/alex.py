import json

def sum_second_elements(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)

    return [sum(item[1] for item in sublist) for sublist in data['coefficient_q_powers']]
