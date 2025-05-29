import sys
import re
from collections import defaultdict

MIN_IDENTITY = 90.0
MIN_COVERAGE = 0.7
MAX_HITS = 100
DROP_THRESHOLD = 3.0
MIN_HITS_FOR_DROP = 10
MIN_TOP_IDENT = 97
CONSENSUS_THRESHOLD = 0.9

TAXONOMIC_ORDER = "dkpcofgs"

def extract_taxonomy(sseqid):
    tax_match = re.search(r";tax=([a-zA-Z0-9:,._\-]+)", sseqid)
    if not tax_match:
        return []

    raw_tax = tax_match.group(1).split(",")
    cleaned_tax = []

    for entry in raw_tax:
        if ":" not in entry:
            continue
        rank, name = entry.split(":", 1)
        if re.search(r"(_g_sp|_sp|_g)$", name):
            continue
        cleaned_tax.append(f"{rank}:{name}")
    
    return cleaned_tax

def lowest_common_ancestor(taxa_list, consensus_threshold=0.9):
    if not taxa_list:
        return []
    
    lca = []
    num_taxa = len(taxa_list)
    min_len = min(len(t) for t in taxa_list)

    for i in range(min_len):
        level_counts = defaultdict(int)
        for tax in taxa_list:
            level_counts[tax[i]] += 1

        most_common_tax, count = max(level_counts.items(), key=lambda x: x[1])
        if count / num_taxa >= consensus_threshold:
            lca.append(most_common_tax)
        else:
            break

    return lca

def rank_summary(taxa_lines):
    level_counts = defaultdict(int)
    level_order = []

    # First pass to get all used levels
    for line in taxa_lines:
        found_levels = re.findall(r'([a-z]):', line)
        for level in found_levels:
            if level not in level_order:
                level_order.append(level)

    level_order = sorted(set(level_order), key=lambda x: TAXONOMIC_ORDER.index(x) if x in TAXONOMIC_ORDER else len(TAXONOMIC_ORDER))

    for line in taxa_lines:
        found_levels = re.findall(r'([a-z]):', line)
        for level in level_order:
            if level in found_levels:
                level_counts[level] += 1

    print("\nClassifiable counts per taxonomic level:")
    for level in level_order:
        print(f"{level}: {level_counts[level]}")

def main(blast_file, output_file):
    hits = defaultdict(list)
    identities = defaultdict(list)
    total_count = 0
    total_seq = 0

    with open(blast_file) as f:
        for line in f:
            total_count += 1
            fields = line.strip().split("\t")
            if len(fields) < 13:
                continue

            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            align_len = int(fields[3])
            qstart, qend = int(fields[6]), int(fields[7])

            qlen = abs(qend - qstart) + 1
            coverage = align_len / qlen if qlen > 0 else 0.0

            if pident >= MIN_IDENTITY and coverage >= MIN_COVERAGE:
                taxonomy = extract_taxonomy(sseqid)
                if taxonomy:
                    hits[qseqid].append(taxonomy)
                    identities[qseqid].append(pident)

    truncated_hits = {}
    for qid, tax_list in hits.items():
        score_list = identities[qid]
        cut_index = len(score_list)

        for i in range(1, min(len(score_list), MAX_HITS)):
            if i >= MIN_HITS_FOR_DROP:
                drop = score_list[i - 1] - score_list[i]
                if drop >= DROP_THRESHOLD:
                    cut_index = i
                    break
        truncated_hits[qid] = tax_list[:cut_index]

    output_lines = []

    with open(output_file, "w") as out:
        for qid in sorted(truncated_hits):
            total_seq += 1

            subset = truncated_hits[qid]
            if not subset:
                out.write(f"{qid}\t\n")
                continue

            lca = lowest_common_ancestor(subset)
            top_identity = identities[qid][0] if identities[qid] else 0

            if top_identity < MIN_TOP_IDENT:
                lca = [tax for tax in lca if not tax.startswith("s:")]

            tax_string = ",".join(lca)
            out.write(f"{qid}\t{tax_string}\n")
            output_lines.append(tax_string)

    # Show rank summary
    print(f"total IDs: {total_count}")
    print(f"total Seqs: {total_seq}")
    rank_summary(output_lines)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python infer_lca.py <blast_output.tsv> <output.tsv>")
        sys.exit(1)

    main(*sys.argv[1:])
