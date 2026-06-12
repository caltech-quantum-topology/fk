from fkcompute.inversion.api import find_sign_assignment_full
import yaml
from pathlib import Path
import copy

braids = [
    [-1,-2,-2,-2,-1,2,2,2]
]

weights = [None]

sign_assignments = []
for braid,weight in zip(braids,weights):
    signs = find_sign_assignment_full(braid, weight=weight)
    for sign in signs:
        sign_assignments.append({"braid": copy.deepcopy(braid), "inversion": sign.sign_assignment})

output_file = Path("output.yaml")
with output_file.open(mode = "w") as f:
    yaml.safe_dump(sign_assignments, f)
