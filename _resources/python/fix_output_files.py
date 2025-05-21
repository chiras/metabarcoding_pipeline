import argparse
import re
import fileinput
import shutil

# Define patterns and replacements as a dictionary
patterns_asv = {
    'pattern1': (r'^#OTU ID', '')
}
patterns_taxonomy = {
    'pattern1': (r',([A-Za-z0-9_-]*;tax=|	)', ','),
    'pattern1.5': (r'\t', ','),
    'pattern2': (r';size=[0-9]*', ''),
    'pattern3': (r',.*d:', ',d:'),
   # 'pattern4': (r'c:.*,o:', 'o:'),
    'pattern5': (r'_[0-9]*,', ','),
    'pattern6': (r',.*k:', ',k:'),
    'pattern7': (r';$', '')
}

def process_file(filename, patterns_dict):
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            modified_line = line
            for pattern, replacement in patterns_dict.values():
                modified_line = re.sub(pattern, replacement, modified_line)
            print(modified_line, end='')

def main():
    parser = argparse.ArgumentParser(description='Script to modify files.')
    parser.add_argument('--asv', help='Specify the path for asv_table.merge.txt')
    parser.add_argument('--tax', help='Specify the path for taxonomy.vsearch')
    
    args = parser.parse_args()
    
    if args.asv:
        process_file(args.asv, patterns_asv)
    
    if args.tax:
        process_file(args.tax, patterns_taxonomy)

if __name__ == "__main__":
    main()
