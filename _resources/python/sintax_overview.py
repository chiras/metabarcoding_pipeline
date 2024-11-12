import re
import sys
from collections import defaultdict

def count_classifiable_levels(filename):
    counts = defaultdict(int)
    unique_levels = set()

    # First pass to identify all unique levels in the file
    with open(filename, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue

            classifications = fields[-1]  # Taking the last column as the classification part
            found_levels = re.findall(r'([a-z]):[A-Za-z_]+\(.*?\)', classifications)
            unique_levels.update(found_levels)
    
    # Sort the levels in a standard taxonomic order, if possible
    ordered_levels = sorted(unique_levels, key=lambda x: "dkgpcofgs".index(x) if x in "dkgpcofgs" else len("dkgpcofgs"))

    # Second pass to count sequences classifiable to each level
    with open(filename, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue

            classifications = fields[-1]  # Taking the last column as the classification part
            found_levels = re.findall(r'([a-z]):[A-Za-z_]+\(.*?\)', classifications)
            for level in ordered_levels:
                if level in found_levels:
                    counts[level] += 1

    # Display results
    print("Classifiable counts per taxonomic level:")
    for level in ordered_levels:
        print(f"{level}: {counts[level]}")


if len(sys.argv) != 2:
    print("Usage: python sintax_overview.py <filename>")
    sys.exit(1)

filename = sys.argv[1]
count_classifiable_levels(filename)
