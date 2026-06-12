import csv
import gzip

def find_where(lst,predicate):
    for x in lst:
        if predicate(x):
            return x
        
def find_index(lst,predicate):
    for i,x in enumerate(lst):
        if predicate(x):
            return i

def csv_to_dicts(csv_file_path):
    data_list = []
    
    # Open file as either gzipped or plain text
    open_func = gzip.open if csv_file_path.endswith('.gz') else open
    
    with open_func(csv_file_path, 'rt') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            data_list.append(row)
    
    return data_list



def tsv_to_dicts(csv_file_path):
    data_list = []
    with open(csv_file_path, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter='\t')
        for row in csv_reader:
            data_list.append(row)
    return data_list

def save_dicts_to_tsv(data, file_path):
    with open(file_path, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=data[0].keys(), delimiter='\t')
        writer.writeheader()
        for row in data:
            writer.writerow(row)

def sort_any(xs):
    return list(sorted(xs, key=lambda x:str(x)))
