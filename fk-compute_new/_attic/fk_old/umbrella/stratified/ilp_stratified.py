import os
import json
from braidstates_stratified import BraidStates
from fkstratified import FK

# regarding user features: ilp out files can be stored with degree specified as free variable, for reuse with different degrees without recomputation

DEGREE = 10

braids = [
    [1, 1, 2, -1, -3, 2, -3, -4, 3, -4],
    [1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 2, 2, 1 ,-2],
    [1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 2, 2, 1, -2],
    [1, 1, 1, 2, 2, 1, -2, 2, 1, -2],
    [1, 3, 2, -3, 3, 2, -3, 3, 2, 1, -2, -3, 3, 2, 1, -2, -3, 2, 3]
]
titles = ["6_1", "3_1", "5_1", "5_2", "7_1", "7_3", "7_5", "8_15"]

import os
for item in list(zip(titles, braids)):
    FK(
        braid=item[1],
        degree=DEGREE,
        weight=5,
        inversion_data={0: [1 for _ in range(2 * len(item[1]))]},
        outfile=f'Data/Input/{item[0]}s.csv'
    )
    os.system(f'./stratified "Data/Input/{item[0]}s" "Data/Output/{item[0]}s"')
