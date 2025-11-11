import requests
import json

url = "https://topology.fyi/api/fk_table_homogeneous"

headers = {
    "Accept": "application/json",
    "Authorization": "Basic dG9wb2xvZ3k6Znlp"
}

print(f'loading data from {url}...')

response = requests.get(url, headers=headers)
data = response.json()
with open("homo_fk.json","w") as f:
    json.dump(data,f)

print(f'...downloaded {len(data)} rows')
