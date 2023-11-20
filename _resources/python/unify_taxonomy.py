import argparse
import re
import os
from Bio import SeqIO

def read_reference_taxonomy(ref_file,target_ranks,synonym_dict):
    ref_taxonomy_dict = {}
    with open(ref_file, 'r') as ref_fasta:
        for line in ref_fasta:
            if line.startswith('>'):
                header = line[1:].strip()
                id, taxonomy = header.split(';tax=')
                # Remove trailing semicolon
                taxonomy = taxonomy.rstrip(';')
                # Replace non-alphabetic characters with a placeholder
                taxonomy = re.sub(r'[^a-zA-Z-_,:]', '-', taxonomy)
                
                for syn, target in synonym_dict.items():
                    taxonomy = taxonomy.replace(syn + ":", target + ":")

                # Split taxonomy into ranks
                ranks = {rank.split(':')[0]: rank.split(':')[1].strip('-') for rank in taxonomy.split(',') if ':' in rank}

                if 'g' in ranks and 'Candidatus' in ranks['g']:
                    ranks['g'] = ranks['g'].replace('_', '-')

                if 'g' in ranks and '_' in ranks['g']:
                    ranks['g'] = ranks['g'].split('_')[0]

                modified_taxonomy = f'' + ','.join([f'{key}:{value}' for key, value in ranks.items() if key in target_ranks])
                # Store g rank and full taxonomy in the dictionary if not already present
                if 'g' in ranks and set(target_ranks).issubset(ranks) and ranks['g'] not in ref_taxonomy_dict:
                    ref_taxonomy_dict[ranks['g']] = modified_taxonomy

    return ref_taxonomy_dict

def read_name_standardizations(file_path):
    name_std_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            bad_name, replacement = map(str.strip, line.split(';'))
            name_std_dict[bad_name] = replacement
    return name_std_dict


def unify_taxonomy(input_file, target_ranks, synonyms, ref_taxonomy_file, name_standardizations):
    synonym_dict = {syn[0]: syn[1] for syn in synonyms}
    ref_taxonomy_dict = read_reference_taxonomy(ref_taxonomy_file, target_ranks,synonym_dict) if ref_taxonomy_file else {}
    name_std_dict = read_name_standardizations(name_standardizations) if name_standardizations else {}

    genus_taxonomy_dict = {}

    output_file_path = input_file.rsplit('.', 1)[0] + '.ut.fa'

    # Remove the existing output file if it exists
    if os.path.exists(output_file_path):
        os.remove(output_file_path)

    with open(input_file, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.id
            id, taxonomy = header.split(';tax=')

            # Remove trailing semicolon
            taxonomy = taxonomy.rstrip(';')

            # Replace non-alphabetic characters with a placeholder
            taxonomy = re.sub(r'[^a-zA-Z-_,:]', '-', taxonomy)

            # Replace synonyms
            for syn, target in synonym_dict.items():
                taxonomy = taxonomy.replace(syn + ":", target + ":")

            # Check if name_std_dict is not empty before attempting to replace names
            if name_std_dict:
                for bad_name, replacement in name_std_dict.items():
                    taxonomy = taxonomy.replace(bad_name, replacement)

            # Split taxonomy into ranks
            ranks = {rank.split(':')[0]: rank.split(':')[1].strip('-') for rank in taxonomy.split(',') if ':' in rank}

            # Check if 'g' is present and all target ranks are present in the taxonomy
            if 'g' in ranks and set(target_ranks).issubset(ranks):

                # Remove everything after "_" in the genus

                if 'g' in ranks and 'Candidatus' in ranks['g']:
                    ranks['g'] = ranks['g'].replace('_', '-')

                if 'g' in ranks and '_' in ranks['g']:
                    ranks['g'] = ranks['g'].split('_')[0]

                # Remove the record if any rank value is "unidentified"
                if any(any(specified.lower() in value.lower() for specified in {'unidentified', 'clone', 'isolate', 'culture', 'strain', 'environmental', 'symbiont'}) for value in ranks.values()):
                    continue

                # Remove "_-" in the family
                if 'f' in ranks and '_-' in ranks['f']:
                    ranks['f'] = ranks['f'].replace('_-', '').replace('_', '')

                # Use ref_taxonomy for g rank if available, otherwise use modified_header
                modified_header = f'>{id};tax=' + ','.join([f'{key}:{value}' for key, value in ranks.items() if key in target_ranks])

                ref_taxonomy = ref_taxonomy_dict.get(ranks['g'], None)

                if ref_taxonomy:
                    # if "Planktothrix" in ref_taxonomy: 
                    #     print(f"ref hit: {modified_header} -- {ref_taxonomy}")
                    # Replace g rank with ref_taxonomy
                    modified_header = f'>{id};tax={ref_taxonomy}'
                else:
                    if ranks['g'] not in genus_taxonomy_dict:
                        genus_taxonomy_dict[ranks['g']] = modified_header
                    else:
                        modified_header = genus_taxonomy_dict[ranks['g']] 

                # if "Planktothrix" in modified_header: 
                #     print(f"ref hit: {modified_header} -- {ref_taxonomy}")

                # Replace non-ACGT characters in the sequence with 'N'
                modified_sequence = re.sub(r'[^ACGTacgt]', 'N', str(record.seq))

                # Remove leading and trailing 'N's from the sequence
                modified_sequence = modified_sequence.lstrip('N').rstrip('N')

                # Save modified header and sequence to a new file
                with open(output_file_path, 'a') as output_file:
                    output_file.write(modified_header + '\n')
                    output_file.write(modified_sequence + '\n')

def main():
    parser = argparse.ArgumentParser(description='Unify Taxonomy in FASTA Headers')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-t', '--ranks', type=str, required=True, help='Target ranks as a string, e.g., "dpcofg"')
    parser.add_argument('-m', '--synonyms', action='append', nargs=2, metavar=('SYN', 'TARGET'),
                        help='Specify synonym pairs, e.g., "-m d k" (for domain -> kingdom))')
    parser.add_argument('--ref', type=str, help='Reference taxonomy file (optional)')
    parser.add_argument('-n', '--name-standardizations-file', type=str, required=False,
                        help='File containing bad names and their replacements in the format BADNAME;REPLACEMENT')
    args = parser.parse_args()

    unify_taxonomy(args.input, args.ranks, args.synonyms, args.ref,args.name_standardizations_file)

if __name__ == "__main__":
    main()

# create standardization file
# ❯ grep ">" ./_DBs/16S_rdp_16s_v18.ut.fa > ./_DBs/16S_headers.txt
# ❯ grep ">" ./_DBs/16S_silva_16s_v123.fa >> ./_DBs/16S_headers.txt
# ❯ grep ">" ./_DBs/16S_gg_16s_13.5.fa >> ./_DBs/16S_headers.txt