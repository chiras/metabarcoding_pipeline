#!/usr/bin/env python3

import argparse
import csv
import os
import re
import subprocess
from collections import defaultdict


def safe_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", x)


def strip_quotes(x):
    x = x.strip()
    if len(x) >= 2 and ((x[0] == '"' and x[-1] == '"') or (x[0] == "'" and x[-1] == "'")):
        return x[1:-1]
    return x


def remove_inline_comment(line):
    """
    Remove simple bash inline comments outside quotes.
    This is intentionally conservative and sufficient for config lines like:
      refDB_path[0]="/path/file.fa"   # comment
    """
    in_single = False
    in_double = False

    for i, ch in enumerate(line):
        if ch == "'" and not in_double:
            in_single = not in_single
        elif ch == '"' and not in_single:
            in_double = not in_double
        elif ch == "#" and not in_single and not in_double:
            return line[:i].rstrip()

    return line.rstrip()


def parse_bash_config_arrays(config_path):
    """
    Parser for bash-3-compatible mixed-marker config style:

      marker_groups[0]="fITS"
      marker_groups[1]="16S"

      refDB_marker[0]=0
      refDB_path[0]="/path/to/fITS_db.fa"

      refDB_marker[1]=1
      refDB_path[1]="/path/to/16S_db1.fa"

    It also tolerates old direct DB style for single-marker configs:

      marker="fITS"
      refDBs[0]="/path/to/db.fa"

    Return:
      markers: list of tuples [(marker_index, marker_name), ...]
      refdbs:  dict marker_index -> dict db_order -> db_path
    """

    marker_groups = {}
    refdb_marker = {}
    refdb_path = {}

    old_marker = None
    old_refdbs = {}

    re_marker_group = re.compile(r'^\s*marker_groups\[(\d+)\]\s*=\s*(.+?)\s*$')
    re_ref_marker = re.compile(r'^\s*refDB_marker\[(\d+)\]\s*=\s*(\d+)\s*$')
    re_ref_path = re.compile(r'^\s*refDB_path\[(\d+)\]\s*=\s*(.+?)\s*$')

    # Backward-compatible old-style single-marker parsing
    re_old_marker = re.compile(r'^\s*marker\s*=\s*(.+?)\s*$')
    re_old_ref = re.compile(r'^\s*refDBs\[(\d+)\]\s*=\s*(.+?)\s*$')

    with open(config_path) as f:
        for raw_line in f:
            line = remove_inline_comment(raw_line).strip()

            if not line:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("declare "):
                continue

            m = re_marker_group.match(line)
            if m:
                marker_groups[int(m.group(1))] = strip_quotes(m.group(2))
                continue

            m = re_ref_marker.match(line)
            if m:
                refdb_marker[int(m.group(1))] = int(m.group(2))
                continue

            m = re_ref_path.match(line)
            if m:
                refdb_path[int(m.group(1))] = strip_quotes(m.group(2))
                continue

            m = re_old_marker.match(line)
            if m:
                old_marker = strip_quotes(m.group(1))
                continue

            m = re_old_ref.match(line)
            if m:
                old_refdbs[int(m.group(1))] = strip_quotes(m.group(2))
                continue

    # New mixed-marker / generic style
    if marker_groups:
        refdbs = defaultdict(dict)

        all_ref_indices = sorted(set(refdb_marker) | set(refdb_path))

        for ref_i in all_ref_indices:
            if ref_i not in refdb_marker:
                raise SystemExit(f"ERROR: refDB_marker[{ref_i}] missing in config")
            if ref_i not in refdb_path:
                raise SystemExit(f"ERROR: refDB_path[{ref_i}] missing in config")

            marker_i = refdb_marker[ref_i]
            if marker_i not in marker_groups:
                raise SystemExit(
                    f"ERROR: refDB_marker[{ref_i}]={marker_i}, "
                    f"but marker_groups[{marker_i}] is not defined"
                )

            # Preserve priority order within each marker by insertion index.
            # The ref_i itself is used as db order. This is OK because sorting by
            # ref_i preserves the order in the config.
            refdbs[marker_i][ref_i] = refdb_path[ref_i]

        markers = [(i, marker_groups[i]) for i in sorted(marker_groups)]

        return markers, refdbs

    # Old single-marker fallback
    if old_refdbs:
        marker_name = old_marker or "marker"
        markers = [(0, marker_name)]
        refdbs = defaultdict(dict)

        for db_i in sorted(old_refdbs):
            refdbs[0][db_i] = old_refdbs[db_i]

        return markers, refdbs

    raise SystemExit(
        "ERROR: no marker_groups[] / refDB_marker[] / refDB_path[] found, "
        "and no old-style refDBs[] fallback found in config"
    )


def fasta_ids(fasta):
    ids = []

    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0].split(";")[0])

    return ids


def run_vsearch(vsearch, query, db, out_uc, threshold, threads, log):
    cmd = [
        vsearch,
        "--usearch_global", query,
        "--db", db,
        "--id", f"0.{threshold}",
        "--uc", out_uc,
        "--strand", "both",
        "--threads", str(threads),
    ]

    with open(log, "w") as err:
        subprocess.run(cmd, check=True, stderr=err)


def parse_uc_hits(uc_file, marker_name, db_path, db_index):
    hits = []

    with open(uc_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if parts[0] != "H":
                continue

            identity = float(parts[3])
            qid = parts[8].split()[0].split(";")[0]
            target = parts[9] if len(parts) > 9 else ""

            hits.append({
                "asv_id": qid,
                "marker": marker_name,
                "identity": identity,
                "db_index": db_index,
                "db": db_path,
                "target": target,
            })

    return hits


def main():
    ap = argparse.ArgumentParser(
        description="Assign ASVs to marker groups by searching ASVs against marker-specific reference DB groups."
    )
    ap.add_argument("--asv-fasta", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--vsearch", required=True)
    ap.add_argument("--threshold", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--outdir", default="marker_split")
    ap.add_argument("--min-delta", type=float, default=2.0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs("logs/marker_split", exist_ok=True)

    markers, refdbs = parse_bash_config_arrays(args.config)

    if len(markers) < 1:
        raise SystemExit("ERROR: no marker groups found in config")

    for marker_i, marker_name in markers:
        if marker_i not in refdbs or not refdbs[marker_i]:
            raise SystemExit(
                f"ERROR: marker group {marker_i} ({marker_name}) has no refDB_path entries"
            )

    all_asvs = fasta_ids(args.asv_fasta)
    all_hits = []

    print("-- marker groups detected")
    for marker_i, marker_name in markers:
        print(f"marker_group[{marker_i}] = {marker_name}")
        for db_i in sorted(refdbs.get(marker_i, {})):
            print(f"  refDB {db_i}: {refdbs[marker_i][db_i]}")

    for marker_i, marker_name in markers:
        marker_safe = safe_name(marker_name)

        for db_i in sorted(refdbs.get(marker_i, {})):
            db = refdbs[marker_i][db_i]

            uc = os.path.join(args.outdir, f"{marker_safe}.db{db_i}.uc")
            log = os.path.join("logs", "marker_split", f"{marker_safe}.db{db_i}.log")

            print(f"-- marker assignment search: marker={marker_name} db={db_i}")
            print(f"DB: {db}")

            run_vsearch(
                args.vsearch,
                args.asv_fasta,
                db,
                uc,
                args.threshold,
                args.threads,
                log,
            )

            all_hits.extend(parse_uc_hits(uc, marker_name, db, db_i))

    hits_by_asv = defaultdict(list)

    all_hits_path = os.path.join(args.outdir, "all_marker_hits.tsv")
    with open(all_hits_path, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=["asv_id", "marker", "identity", "db_index", "db", "target"],
            delimiter="\t",
        )
        writer.writeheader()

        for hit in all_hits:
            writer.writerow(hit)
            hits_by_asv[hit["asv_id"]].append(hit)

    assigned = defaultdict(list)

    assignment_path = "marker_assignment.tsv"
    assignment_path_in_outdir = os.path.join(args.outdir, "marker_assignment.tsv")

    with open(assignment_path, "w", newline="") as out, \
         open(assignment_path_in_outdir, "w", newline="") as out2:

        fieldnames = [
            "asv_id",
            "assigned_marker",
            "best_identity",
            "second_identity",
            "delta",
            "best_db_index",
            "best_db",
            "status",
        ]

        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer2 = csv.DictWriter(out2, fieldnames=fieldnames, delimiter="\t")

        writer.writeheader()
        writer2.writeheader()

        for asv in all_asvs:
            if asv not in hits_by_asv:
                assigned["unassigned"].append(asv)

                row = {
                    "asv_id": asv,
                    "assigned_marker": "unassigned",
                    "best_identity": "NA",
                    "second_identity": "NA",
                    "delta": "NA",
                    "best_db_index": "NA",
                    "best_db": "NA",
                    "status": "unassigned",
                }

                writer.writerow(row)
                writer2.writerow(row)
                continue

            best_by_marker = {}

            for hit in hits_by_asv[asv]:
                marker = hit["marker"]

                if marker not in best_by_marker:
                    best_by_marker[marker] = hit
                elif hit["identity"] > best_by_marker[marker]["identity"]:
                    best_by_marker[marker] = hit

            ranked = sorted(
                best_by_marker.values(),
                key=lambda x: x["identity"],
                reverse=True,
            )

            best = ranked[0]
            second_identity = "NA"
            delta = "NA"

            if len(ranked) == 1:
                assigned_marker = best["marker"]
                status = "assigned"
            else:
                second_identity = ranked[1]["identity"]
                delta = best["identity"] - ranked[1]["identity"]

                if delta < args.min_delta:
                    assigned_marker = "ambiguous"
                    status = "ambiguous"
                else:
                    assigned_marker = best["marker"]
                    status = "assigned"

            assigned[assigned_marker].append(asv)

            row = {
                "asv_id": asv,
                "assigned_marker": assigned_marker,
                "best_identity": best["identity"],
                "second_identity": second_identity,
                "delta": delta,
                "best_db_index": best["db_index"],
                "best_db": best["db"],
                "status": status,
            }

            writer.writerow(row)
            writer2.writerow(row)

    # Write ID files for all marker groups, including empty files.
    for _, marker_name in markers:
        marker_safe = safe_name(marker_name)

        with open(os.path.join(args.outdir, f"{marker_safe}.ids"), "w") as out:
            for asv in assigned.get(marker_name, []):
                out.write(asv + "\n")

    for special in ["ambiguous", "unassigned"]:
        with open(os.path.join(args.outdir, f"{special}.ids"), "w") as out:
            for asv in assigned.get(special, []):
                out.write(asv + "\n")

    print("-- marker assignment summary")
    total_assigned = 0

    for _, marker_name in markers:
        n = len(assigned.get(marker_name, []))
        total_assigned += n
        print(f"{marker_name}: {n}")

    n_ambiguous = len(assigned.get("ambiguous", []))
    n_unassigned = len(assigned.get("unassigned", []))

    print(f"ambiguous: {n_ambiguous}")
    print(f"unassigned: {n_unassigned}")
    print(f"total ASVs: {len(all_asvs)}")
    print(f"assigned to marker groups: {total_assigned}")


if __name__ == "__main__":
    main()