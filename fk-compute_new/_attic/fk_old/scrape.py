import requests
import json
import os

url = "https://topology.fyi/api/Verified%20sign%20diagrams"

headers = {
    "Accept": "application/json",
    "Authorization": "Basic dG9wb2xvZ3k6Znlp"
}

print(f'loading data from {url}...')

response = requests.get(url, headers=headers)
data = response.json()

print(f'...downloaded {len(data)} rows')

# Create the directory if it doesn't exist
os.makedirs('Data/Input', exist_ok=True)

# Write the data to the JSON file
with open('Data/Input/inversion_data.json', 'w') as f:
    json.dump(data, f, indent=2)

print('Data saved to Data/Input/inversion_data.json')
