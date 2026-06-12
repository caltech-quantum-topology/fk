import csv

def string_to_list(s):
    return [int(x) for x in s.strip('{}').split(';')]

with open('/home/josefsvoboda/Mathematics/projects/quantum_invariants_and_fibered_knots/code/fast_fk_code/fk1/fk/2bridge_knots_9_cross.csv', newline='') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    titres = []
    braids = []
    next(data)
    for row in data:
       titres.append(row[0])
       braids.append(string_to_list(row[5]))
    print(braids)
