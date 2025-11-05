#!/bin/bash

# Count total lines of .cpp, .hpp, and Makefile files recursively

total=0

# Use find to locate all matching files (case-insensitive for makefile)
while IFS= read -r file; do
    lines=$(wc -l < "$file")
    total=$((total + lines))
done < <(find . -type f \( -iname "*.cpp" -o -iname "*.hpp" -o -iname "Makefile" \))

echo "Total lines of code (cpp/hpp/Makefiles): $total"
