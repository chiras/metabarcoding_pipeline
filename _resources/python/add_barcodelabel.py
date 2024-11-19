import re
import argparse

parser = argparse.ArgumentParser(description="Process FASTA file and add barcodelabels.")
parser.add_argument("-i", "--input", required=True, help="Path to the input FASTA file")
parser.add_argument("-o", "--output", required=True, help="Path to the output FASTA file")
args = parser.parse_args()

with open(args.input, "r") as infile, open(args.output, "w") as outfile:
    for line in infile:
        if line.startswith(">R1+2-"):
            line = re.sub(r"^>R1\+2-(.*)_(\d+);", r">R1+2-\1_\2;barcodelabel=\1;", line)
        elif line.startswith(">R1-"):
            line = re.sub(r"^>R1-(.*)_(\d+)$", r">R1-\1_\2;barcodelabel=\1", line)
        outfile.write(line)
