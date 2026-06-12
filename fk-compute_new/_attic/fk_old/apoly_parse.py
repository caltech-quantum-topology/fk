import os
import re

# Setup
input_dir = "data_apoly"
output_file = "polynomials.m"

# Get files starting with L108, L109, or L110
all_files = [f for f in os.listdir(input_dir) if f.endswith('.txt')]
selected_files = sorted([f for f in all_files 
                  if f.startswith('L108') or f.startswith('L109') or f.startswith('L110')])

print(f"Found {len(selected_files)} files to process")

# Start writing Mathematica file
with open(output_file, 'w') as out:
    out.write("<|\n")  # Start Association
    
    for i, filename in enumerate(selected_files, 1):
        filepath = os.path.join(input_dir, filename)
        
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Split on :=
        if ':=' not in content:
            continue
            
        parts = content.split(':=', 1)
        var_name = parts[0].strip()
        expression = parts[1].strip()
        
        # Remove semicolons
        expression = expression.replace(';', '')
        
        # Remove line breaks and extra whitespace
        expression = expression.replace('\n', ' ').replace('\r', ' ')
        expression = re.sub(r'\s+', ' ', expression)
        expression = expression.strip()
        
        # Write as Mathematica association entry: "key" -> value,
        out.write(f'  "{var_name}" -> {expression}')
        
        # Add comma if not last entry
        if i < len(selected_files):
            out.write(',\n')
        else:
            out.write('\n')
        
        if i % 100 == 0:
            print(f"Processed {i}/{len(selected_files)}")
    
    out.write("|>")  # End Association

print(f"Done! Mathematica code written to {output_file}")
