import argparse
import re
from collections import defaultdict
import os

if not os.path.exists("logs"):
    os.makedirs("logs")

TAXONOMIC_RANKS = ["d", "k", "p", "c", "o", "f", "g", "s"]

def parse_taxonomy_string(tax_string):
    """Convert taxonomy string into dict of rank:taxon"""
    tax_dict = {}
    for part in tax_string.strip().split(","):
        if ":" in part:
            k, v = part.strip().split(":", 1)
            tax_dict[k] = v
    return tax_dict

def detect_rank_order(tax_lines):
    """Detect rank order from a list of taxonomy lines"""
    for line in tax_lines:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            tax = parts[-1]
            ranks = []
            for part in tax.strip().split(","):
                match = re.match(r"([a-z]):[^,:]+", part)
                if match:
                    rank = match.group(1)
                    if rank not in ranks:
                        ranks.append(rank)
            if ranks:
                return ranks
    return []  # fallback

def merge_taxonomies(asv, blast_tax, sintax_tax, w_blast=0.6, w_sintax=0.4, greedy=False):    
    final_tax = {}
    all_ranks = [r for r in TAXONOMIC_RANKS if r in blast_tax or r in sintax_tax]

    for r in all_ranks:
        b_val = blast_tax.get(r)
        s_val = sintax_tax.get(r)

        if b_val == s_val and b_val is not None:
            final_tax[r] = b_val
        elif greedy:
            # Take whichever is filled if only one is filled
            if b_val and not s_val:
                final_tax[r] = b_val
            elif s_val and not b_val:
                final_tax[r] = s_val
            # If both filled and disagree, skip this rank
        elif b_val and s_val:
            final_tax[r] = b_val if w_blast >= w_sintax else s_val
        elif b_val:
            final_tax[r] = b_val
        elif s_val:
            final_tax[r] = s_val

    return final_tax



def format_tax_dict(tax_dict, rank_order):
    """Format taxonomy dictionary using detected rank order"""
    return ",".join(f"{r}:{tax_dict[r]}" for r in rank_order if r in tax_dict)

def main(blast_file, sintax_file, output_file, w_blast, w_sintax, greedy):    
    blast = {}
    sintax = {}
    all_lines = []

    # Read BLAST taxonomy
    with open(blast_file) as f:
        blast_lines = f.readlines()
        for line in blast_lines:
            if not line.strip():
                continue
            parts = line.strip().split("\t", 1)
            if len(parts) < 2:
                continue
            asv, tax = parts
            blast[asv] = parse_taxonomy_string(tax)

    # Read SINTAX taxonomy
    with open(sintax_file) as f:
        sintax_lines = f.readlines()
        for line in sintax_lines:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            asv = fields[0]
            tax = fields[-1]
            if not tax.strip():
                continue
            sintax[asv] = parse_taxonomy_string(tax)

    # Detect rank order from combined sources
    #rank_order = detect_rank_order(blast_lines + sintax_lines)
    rank_order = TAXONOMIC_RANKS
    # Merge
    all_asvs = set(blast.keys()) | set(sintax.keys())
    unmatched = []
    detailed_report = []

    with open(output_file, "w") as out:
        for asv in sorted(all_asvs):
            b_tax = blast.get(asv, {})
            s_tax = sintax.get(asv, {})
            if not b_tax and not s_tax:
                unmatched.append(asv)
                continue
            merged = merge_taxonomies(asv, b_tax, s_tax, w_blast, w_sintax, greedy)
            final_str = format_tax_dict(merged, rank_order)
            # Replace ';size=digits' in ASV ID + TAB with comma
            # Example: ASV10000;size=258 \t k:Animalia,...
            asv_clean = re.sub(r";size=\d+", "", asv)
            out.write(f"{asv_clean},{final_str}\n")


    # Report
    print("=== Merge Summary ===")
    print(f"Total ASVs: {len(all_asvs)}")
    print(f"  - Found in BLAST only: {len(set(blast) - set(sintax))}")
    print(f"  - Found in SINTAX only: {len(set(sintax) - set(blast))}")
    print(f"  - Found in both: {len(set(blast) & set(sintax))}")
    print(f"  - Unmatched in both: {len(unmatched)}")
    print()
   # print("=== Sample Merged Entries ===")
  #  for asv, b_len, s_len, m_len in detailed_report[:10]:
   #     print(f"{asv}: BLAST={b_len} ranks, SINTAX={s_len} ranks -> Combined={m_len}")

    blast_only = set(blast) - set(sintax)
    sintax_only = set(sintax) - set(blast)
    both = set(blast) & set(sintax)
    #unmatched = set()  # you have an unmatched list from before; keep it as a set for convenience

    with open("logs/blastLCA_sintax.blast_only.txt", "w") as f:
        for asv in sorted(blast_only):
            f.write(asv + "\n")

    with open("logs/blastLCA_sintax.sintax_only.txt", "w") as f:
        for asv in sorted(sintax_only):
            f.write(asv + "\n")

    with open("logs/blastLCA_sintax.both.txt", "w") as f:
        for asv in sorted(both):
            f.write(asv + "\n")

    with open("logs/blastLCA_sintax.unmatched.txt", "w") as f:
        for asv in sorted(unmatched):
            f.write(asv + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", required=True, help="BLAST LCA taxonomy file")
    parser.add_argument("--sintax", required=True, help="SINTAX taxonomy file (cut -f1,4)")
    parser.add_argument("--output", required=True, help="Output merged taxonomy file")
    parser.add_argument("--weight_blast", type=float, default=0.6, help="Weight for BLAST")
    parser.add_argument("--weight_sintax", type=float, default=0.4, help="Weight for SINTAX")
    parser.add_argument("--greedy", action="store_true", help="Use greedy merging (assign taxonomy as deep as possible)")
    args = parser.parse_args()

    main(args.blast, args.sintax, args.output, args.weight_blast, args.weight_sintax, args.greedy)