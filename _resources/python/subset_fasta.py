from tqdm import tqdm
import argparse

def read_text_file(text_file):
    with open(text_file, 'r') as f:
        return set(line.strip() for line in f)

def subset_fasta(input_fasta, text_file, output_fasta):
    header_lines = []
    sequences = []
    with open(input_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header_lines.append(line.strip())
                sequences.append('')
            else:
                sequences[-1] += line.strip()
    
    # Read text from the text file
    target_text = read_text_file(text_file)

    # Filter sequences based on whether the header includes text from the text file
    filtered_headers = []
    filtered_sequences = []
    with tqdm(total=len(header_lines), desc="Processing") as pbar:
        for header, sequence in zip(header_lines, sequences):
            pbar.update(1)
            for text in target_text:
                if text in header:
                    filtered_headers.append(header)
                    filtered_sequences.append(sequence)
                    break

    # Write the subset to the output file
    with open(output_fasta, 'w') as f:
        for header, sequence in zip(filtered_headers, filtered_sequences):
            f.write(header + '\n')
            f.write(sequence + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Subset FASTA file based on text in a text file.')
    parser.add_argument('-i', '--input', dest='input_fasta', required=True, help='Path to the input FASTA file')
    parser.add_argument('-l', '--list', dest='text_file', required=True, help='Path to the text file containing header information')
    parser.add_argument('-o', '--output', dest='output_fasta', required=True, help='Path to the output FASTA file')
    args = parser.parse_args()
    
    subset_fasta(args.input_fasta, args.text_file, args.output_fasta)
