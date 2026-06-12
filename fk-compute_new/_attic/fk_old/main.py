from sign_diagrams import get_sign_assignment, load_sign_data
from braidstates import BraidStates
import alex
import os
import csv
from sympy import symbols, series, collect
import braid_ilp

def load_alexander_vectors(csv_file):
    vectors = {}
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            name = row[0]
            if len(row) > 1 and row[1].strip():  # Check if vector exists and is not empty
                vector_str = row[1].strip('{}')
                if vector_str:  # Additional check for non-empty string
                    vector = [int(x.strip()) for x in vector_str.split(';') if x.strip()]
                    vectors[name] = vector
    return vectors

def compute_reciprocal_polynomial(alex_vector, target_length):
    t = symbols('t')

    # Parse vector format: [min_degree, max_degree, coeff1, coeff2, ...]
    min_degree = alex_vector[0]
    max_degree = alex_vector[1]
    coefficients = alex_vector[2:]

    print(f"Min degree: {min_degree}, Max degree: {max_degree}")
    print(f"Coefficients: {coefficients}")

    # Create polynomial from coefficients starting at min_degree
    poly = sum(coeff * t**(min_degree + i) for i, coeff in enumerate(coefficients))
    print(f"Polynomial: {poly}")

    # Compute reciprocal polynomial: 1/A(t)
    reciprocal_poly = 1/poly
    reciprocal_poly *= -(1 - 1 / t)
    print(f"Reciprocal polynomial: {reciprocal_poly}")

    # First expand the reciprocal polynomial as a series to handle negative powers
    # We need a wide enough range to capture negative powers
    series_expansion = series(reciprocal_poly, t, 0, n=target_length*2).removeO()
    print(f"Series expansion: {series_expansion}")

    # Now collect the expanded series
    collected = collect(series_expansion, t)
    print(f"Collected form: {collected}")

    # Find the minimum degree (could be negative)
    min_power = None
    for i in range(-target_length, target_length):
        coeff = collected.coeff(t, i)
        if coeff is not None and coeff != 0:
            print(f"t^{i}: {coeff}")
            if min_power is None:
                min_power = i

    print(f"Minimum power found: {min_power}")

    # Extract coefficients starting from the minimum power
    coeffs = []
    if min_power is not None:
        for i in range(target_length):
            coeff = collected.coeff(t, min_power + i)
            coeffs.append(int(coeff) if coeff is not None else 0)

    # Remove trailing zeros to match the strategy in compare_vectors
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()

    # Pad or truncate to match target_length
    if len(coeffs) < target_length:
        coeffs.extend([0] * (target_length - len(coeffs)))
    else:
        coeffs = coeffs[:target_length]

    return coeffs

def compare_vectors(knot_id, output_vector, alexander_vectors):
    if knot_id in alexander_vectors:
        alex_vector = alexander_vectors[knot_id]

        # Compute reciprocal polynomial expanded to match output vector length
        reciprocal_coeffs = compute_reciprocal_polynomial(alex_vector, len(output_vector))
        print(reciprocal_coeffs)

        first_nonzero_in_output_intex = next((i for i, x in enumerate(output_vector) if x != 0), None)
        nonzero_output_vector = output_vector[first_nonzero_in_output_intex:]

        print(f"Knot {knot_id}:")
        print(f"  Output vector: {nonzero_output_vector}")
        print(f"  Alexander vector: {alex_vector}")
        print(f"  Reciprocal coeffs: {reciprocal_coeffs}")
        boolean = (nonzero_output_vector == reciprocal_coeffs[:len(nonzero_output_vector)]) or  (list(map(lambda x: -x, nonzero_output_vector)) == reciprocal_coeffs[:len(nonzero_output_vector)])
        print(f"  Match: {boolean}")
        return boolean
    else:
        print(f"Knot {knot_id} not found in Alexander vectors")
        return False
# Load Alexander vectors from CSV
alexander_vectors = load_alexander_vectors('Data/Input/alexander.csv')

def main(title, braid, degree):
    print(f"Processing knot {title} with braid {braid} and degree {degree}")
    assignment = get_sign_assignment(braid, degree=10, verbose=False)
    braid = assignment['braid']
    bs = BraidStates(braid)
    load_sign_data(bs, assignment['inversion_data'], compute_r_matrices=True)
    braid_ilp.save(
        braid_states=bs,
        degree=degree,
        save_to=f'Data/Input/{title}.csv',
    )
    os.system(f'./_ "Data/Input/{title}" "Data/Output/{title}"')

    if not os.path.exists(f'Data/Output/{title}.json'):
        raise Exception(f"Output file for knot {title} not found after C++ computation")

    output_alex_vector = alex.sum_second_elements(f'Data/Output/{title}.json')

    # Compare with Alexander vectors
    if (not compare_vectors(title, output_alex_vector, alexander_vectors)) :
        raise Exception(f"Alexander^-1 mismatch for knot {title}")

if __name__ == "__main__":
    DEGREE = 33
    titres = ["12n_242"]
    braids = [
    [1,2,2,1,1,2,2,2,2,2,2,2],
    ]
    
    for item in list(zip(titres, braids)):
        main(item[0], item[1], DEGREE)
    print("All computations and checks completed successfully.")
