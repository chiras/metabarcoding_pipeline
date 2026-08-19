import argparse
import re
import csv
import fileinput


patterns_asv = {
    "pattern1": (r"^#OTU ID", "")
}


def process_file(filename, patterns_dict):
    with fileinput.FileInput(filename, inplace=True, backup=".bak") as file:
        for line in file:
            modified_line = line
            for pattern, replacement in patterns_dict.values():
                modified_line = re.sub(pattern, replacement, modified_line)
            print(modified_line, end="")


def clean_taxon_name(value):
    """Remove size annotations and trailing numeric ASV artefacts."""
    value = re.sub(r";size=\d+", "", value)
    value = re.sub(r"_\d+$", "", value)
    value = value.strip()
    return value


def parse_taxonomy_string(tax_string):
    """
    Convert a comma- or semicolon-separated taxonomy string into fixed-rank columns.

    Expected rank prefixes:
    d: or k: kingdom/domain
    p: phylum
    c: class
    o: order
    f: family
    g: genus
    s: species
    """

    ranks = {
        "kingdom": "",
        "phylum": "",
        "class": "",
        "order": "",
        "family": "",
        "genus": "",
        "species": ""
    }

    tax_string = tax_string.strip()
    tax_string = re.sub(r";size=\d+", "", tax_string)
    tax_string = tax_string.rstrip(";,")

    # VSEARCH taxonomy may use either commas or semicolons as rank separators.
    parts = [
        p.strip()
        for p in re.split(r"[;,]", tax_string)
        if p.strip()
    ]

    for part in parts:
        part = clean_taxon_name(part)

        if part.startswith(("d:", "k:")):
            ranks["kingdom"] = part
        elif part.startswith("p:"):
            ranks["phylum"] = part
        elif part.startswith("c:"):
            ranks["class"] = part
        elif part.startswith("o:"):
            ranks["order"] = part
        elif part.startswith("f:"):
            ranks["family"] = part
        elif part.startswith("g:"):
            ranks["genus"] = part
        elif part.startswith("s:"):
            ranks["species"] = part

    return ranks


def convert_vsearch_taxonomy_to_csv(infile, outfile):
    """
    Convert VSEARCH taxonomy output to a clean comma-separated table.

    Output columns:
    ASV, kingdom, phylum, class, order, family, genus, species
    """

    fieldnames = [
        "ASV",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
    ]

    with open(infile, "r", encoding="utf-8") as fin, \
         open(outfile, "w", encoding="utf-8", newline="") as fout:

        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()

        for line in fin:
            line = line.strip()

            if not line:
                continue

            # Skip existing taxonomy header.
            if line.startswith(",kingdom,") or line.startswith("ASV,kingdom,"):
                continue

            # Standard VSEARCH output:
            # ASV123;size=100<TAB>k:Viridiplantae,p:...,s:...
            if "\t" in line:
                asv_id, tax_string = line.split("\t", 1)

            # Also support already comma-separated input:
            # ASV123,k:Viridiplantae,p:...,s:...
            elif "," in line:
                asv_id, tax_string = line.split(",", 1)

            else:
                continue

            asv_id = clean_taxon_name(asv_id)
            tax_string = tax_string.strip()

            # Ignore malformed/non-taxonomy lines.
            if not re.search(r"(^|[,;])(d:|k:)", tax_string):
                continue

            ranks = parse_taxonomy_string(tax_string)

            writer.writerow({
                "ASV": asv_id,
                **ranks
            })

def main():
    parser = argparse.ArgumentParser(description="Clean ASV and taxonomy files.")
    parser.add_argument("--asv", help="Path to asv_table.merge.txt")
    parser.add_argument("--tax", help="Path to taxonomy.vsearch")
    parser.add_argument(
        "--tax-out",
        default="taxonomy2.vsearch",
        help="Output path for cleaned taxonomy CSV"
    )

    args = parser.parse_args()

    if args.asv:
        process_file(args.asv, patterns_asv)

    if args.tax:
        convert_vsearch_taxonomy_to_csv(args.tax, args.tax_out)


if __name__ == "__main__":
    main()