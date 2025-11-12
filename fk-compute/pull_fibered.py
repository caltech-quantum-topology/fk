import requests
import json

url = "https://topology.fyi/api/fk_table_fibered"

headers = {
    "Accept": "application/json",
    "Authorization": "Basic dG9wb2xvZ3k6Znlp"
}

print(f'loading data from {url}...')

response = requests.get(url, headers=headers)
data = response.json()
print(data[0].keys())

print(f'...downloaded {len(data)} rows')

filtered_data = []
for row in data:
    filtered_row = {
        "name": row.get("knotinfo_id"),
        "braid": row.get("braid"),
        "max_x_degree": row.get("max_x_degree"),
        "fk_terms": row.get("fk_terms"),
        "fk_x_coeffs": row.get("fk_x_coeffs")
    }
    filtered_data.append(filtered_row)

with open("fibered_table.json", "w") as f:
    json.dump(filtered_data, f, indent=2)

print(f'saved {len(filtered_data)} rows to fibered_table.json')
