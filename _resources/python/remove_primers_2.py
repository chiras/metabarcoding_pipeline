import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import IUPACData
from tqdm import tqdm  # Import tqdm for progress bar

# Define primer presets
PRIMER_PRESETS = {
    "ITS2": {"fw": "ATGCGATACTTGGTGTGAAT", "rv": "TCCTCCGCTTATTGATATGC"},
    "COI": {"fw": "GGWACWGGWTGAACWGTWTAYCCYCC", "rv": "TAAACTTCAGGGTGACCAAAAAATCA"},
    "COI-5P": {"fw": "GGWACWGGWTGAACWGTWTAYCCYCC", "rv": "TAAACTTCAGGGTGACCAAAAAATCA"},
    "16S": {"fw": "CCTACGGGAGGCAGCAG", "rv": "GTGCCAGCMGCCGCGGTAA"}
}

def parse_args():
    parser = argparse.ArgumentParser(description="Remove primers from sequences in a FASTA file.")
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output FASTA file')
    parser.add_argument('--fw', help='Forward primer sequence')
    parser.add_argument('--rv', help='Reverse primer sequence (complement will be removed)')
    parser.add_argument('--marker', choices=PRIMER_PRESETS.keys(), 
                        help='Marker name to use preset forward and reverse primer sequences')
    return parser.parse_args()

def iupac_to_regex(iupac_seq):
    iupac_dict = IUPACData.ambiguous_dna_values
    regex_seq = ''.join([f"[{iupac_dict[nuc]}]" for nuc in iupac_seq])
    return regex_seq

def remove_primers(seq, fw_primer=None, rv_primer=None):
    trimmed = {'start': False, 'end': False}
    seq_str = str(seq)  # Convert to string if it's a Seq object
    
    if fw_primer:
        fw_regex = re.compile(f"^{iupac_to_regex(fw_primer)}", re.IGNORECASE)
        if fw_regex.match(seq_str):
            seq_str = re.sub(fw_regex, "", seq_str)
            trimmed['start'] = True

    if rv_primer:
        rv_seq = str(Seq(rv_primer).reverse_complement())
        rv_regex = re.compile(f"{iupac_to_regex(rv_seq)}$", re.IGNORECASE)
        if rv_regex.search(seq_str):
            seq_str = re.sub(rv_regex, "", seq_str)
            trimmed['end'] = True

    return Seq(seq_str), trimmed

def main():
    args = parse_args()

    # Determine primers from either --fw and --rv or the --marker preset
    if args.marker:
        primers = PRIMER_PRESETS[args.marker]
        fw_primer = primers["fw"]
        rv_primer = primers["rv"]
    elif args.fw or args.rv:
        fw_primer = args.fw
        rv_primer = args.rv
    else:
        raise ValueError("Either both --fw and --rv or a --marker must be provided.")

    trimmed_count = 0
    untrimmed_count = 0

    records = []
    
    # Using tqdm to show a progress bar for processing sequences
    with tqdm(total=sum(1 for _ in SeqIO.parse(args.input, "fasta")), desc="Processing sequences") as pbar:
        for record in SeqIO.parse(args.input, "fasta"):
            new_seq, trimmed = remove_primers(record.seq, fw_primer, rv_primer)
            if trimmed['start'] or trimmed['end']:
                trimmed_count += 1
            else:
                untrimmed_count += 1
            new_record = SeqRecord(new_seq, id=record.id, description=record.description)
            records.append(new_record)
            pbar.update(1)  # Update progress bar for each record

    SeqIO.write(records, args.output, "fasta")

    print(f"Number of sequences trimmed: {trimmed_count}")
    print(f"Number of sequences not trimmed: {untrimmed_count}")

if __name__ == "__main__":
    main()
