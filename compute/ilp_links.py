import os
import json
from braidstates_links import BraidStates
from fklinks import FK





DEGREE = 6





braid = [2, 2, 1, 2, 3, 3, 2, -1, -3, 2, -3, -3]
inversion_data = {0: [-1, 1, 1, 1, -1, -1], 1: [1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1]}

FK(
    braid=braid,
    degree=DEGREE, 
    inversion_data=inversion_data,
    outfile='davide/A.csv'
)

braid = [1, 1, -3, -2, -2, -3, -2, 1, -2, -2]
inversion_data = {0: [1, 1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, -1, -1, 1], 1: [-1, -1, 1, -1]}

FK(
    braid=braid,
    degree=DEGREE, 
    inversion_data=inversion_data,
    outfile='davide/B.csv'
)

braid = [-1, -4, -3, -2, -1, -1, -2, 3, 4, 2, -1, 2, -3, 2, -3]
inversion_data = {0: [-1, -1, -1, -1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1], 1: [-1, -1, -1, -1, -1, -1, -1, -1]}

FK(
    braid=braid,
    degree=DEGREE, 
    inversion_data=inversion_data,
    outfile='davide/C.csv'
)

braid = [-3, -4, -3, -3, 4, 2, -3, 2, -1, 2, -1]
inversion_data = {0: [-1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1], 1: [-1, -1, 1, -1]}

FK(
    braid=braid,
    degree=DEGREE, 
    inversion_data=inversion_data,
    outfile='davide/D.csv'
)

quit()
for file in os.listdir("linksigns"):
    with open(f"linksigns/{file}", 'r') as f:
        data = json.load(f)
    if data["inversion_data"] != "failure":
        linkinfo_id = file[:-5]
        data['inversion_data'] = dict(enumerate(list(data['inversion_data'].values())))
        print(linkinfo_id)
        FK(
            braid=data['braid'],
            degree=4, 
            inversion_data=data['inversion_data'],
            outfile=f'ilp/{linkinfo_id}.csv'
        )