import requests
from bs4 import BeautifulSoup
import os

# Setup
base_url = "https://homepages.math.uic.edu/~culler/Apolynomials/apolys/"
output_dir = "data_apoly"

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Get directory listing
response = requests.get(base_url)
soup = BeautifulSoup(response.text, 'html.parser')

# Find all .apoly files
files = [link.get('href') for link in soup.find_all('a') 
         if link.get('href', '').endswith('.apoly')]

print(f"Found {len(files)} files")

# Download each file
for i, filename in enumerate(files, 1):
    url = base_url + filename
    print(f"Downloading {i}/{len(files)}: {filename}")
    
    # Download content
    content = requests.get(url).text
    
    # Save as .txt file (remove .apoly extension and add .txt)
    txt_filename = filename.replace('.apoly', '.txt')
    filepath = os.path.join(output_dir, txt_filename)
    
    with open(filepath, 'w') as f:
        f.write(content)

print(f"Download complete! All files saved in {output_dir}/")



