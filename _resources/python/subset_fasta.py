#!/usr/bin/env python3

from tqdm import tqdm
import argparse


def normalize_id(value):
    """
    Normalize FASTA/list IDs for exact ASV matching.

    Examples:
      >ASV123;size=45 -> ASV123
      ASV123;size=45  -> ASV123
      ASV123 something -> ASV123
      ASV123          -> ASV123
    """
    value = value.strip()
    if value.startswith(">"):
        value = value[1:]
    value = value.split()[0]
    value = value.split(";")[0]
    return value


def read_text_file(text_file, exact_match=False):
    with open(text_file, "r") as f:
        values = [line.strip() for line in f if line.strip()]

    if exact_match:
        return {normalize_id(value) for value in values}

    # Legacy / greedy behaviour: keep lines as supplied and match by substring.
    return set(values)


def read_fasta(input_fasta):
    header_lines = []
    sequences = []

    with open(input_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                header_lines.append(line.strip())
                sequences.append("")
            else:
                if not sequences:
                    raise ValueError("Input FASTA appears to contain sequence before first header")
                sequences[-1] += line.strip()

    return header_lines, sequences


def subset_fasta(input_fasta, text_file, output_fasta, exact_match=False):
    header_lines, sequences = read_fasta(input_fasta)

    # Read text/IDs from list file.
    target_text = read_text_file(text_file, exact_match=exact_match)

    filtered_headers = []
    filtered_sequences = []

    with tqdm(total=len(header_lines), desc="Processing") as pbar:
        for header, sequence in zip(header_lines, sequences):
            pbar.update(1)

            if exact_match:
                header_id = normalize_id(header)
                if header_id in target_text:
                    filtered_headers.append(header)
                    filtered_sequences.append(sequence)
            else:
                # Legacy / greedy behaviour: any list entry may occur anywhere in header.
                for text in target_text:
                    if text in header:
                        filtered_headers.append(header)
                        filtered_sequences.append(sequence)
                        break

    with open(output_fasta, "w") as f:
        for header, sequence in zip(filtered_headers, filtered_sequences):
            f.write(header + "\n")
            f.write(sequence + "\n")

    print(f"Input FASTA records: {len(header_lines)}")
    print(f"Requested list entries: {len(target_text)}")
    print(f"Written FASTA records: {len(filtered_headers)}")
    print(f"Matching mode: {'exact normalized ID' if exact_match else 'legacy greedy substring'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset FASTA file based on IDs/text in a list file.")
    parser.add_argument("-i", "--input", dest="input_fasta", required=True, help="Path to the input FASTA file")
    parser.add_argument("-l", "--list", dest="text_file", required=True, help="Path to the text file containing IDs/header text")
    parser.add_argument("-o", "--output", dest="output_fasta", required=True, help="Path to the output FASTA file")
    parser.add_argument(
        "--exact_match",
        action="store_true",
        help=(
            "Use exact normalized ID matching instead of legacy greedy substring matching. "
            "Both FASTA headers and list entries are normalized by removing leading '>', "
            "anything after whitespace, and anything after ';'."
        ),
    )

    args = parser.parse_args()

    subset_fasta(
        args.input_fasta,
        args.text_file,
        args.output_fasta,
        exact_match=args.exact_match,
    )
