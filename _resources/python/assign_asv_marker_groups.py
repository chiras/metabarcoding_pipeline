#!/usr/bin/env python3

import argparse
import csv
import os
import re
import subprocess
import shlex
import sys
import time
from collections import defaultdict

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass


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


def parse_bash_config_arrays(config_path, db_source="refDB"):
    """
    Parser for bash-3-compatible mixed-marker config style.

    Supports direct reference DBs and hierarchical/global DBs:

      marker_groups[0]="fITS"
      marker_groups[1]="16S"

      refDB_marker[0]=0
      refDB_path[0]="/path/to/fITS_direct.fa"

      hieDB_marker[0]=0
      hieDB_path[0]="/path/to/fITS_global.fa"

    db_source controls which database set is used for marker assignment:
      refDB  = use refDB_marker/refDB_path
      hieDB  = use hieDB_marker/hieDB_path
      both   = use both sets

    Return:
      markers: list of tuples [(marker_index, marker_name), ...]
      refdbs:  dict marker_index -> dict db_order -> db_path
    """

    marker_groups = {}
    refdb_marker = {}
    refdb_path = {}
    hiedb_marker = {}
    hiedb_path = {}

    old_marker = None
    old_refdbs = {}

    re_marker_group = re.compile(r'^\s*marker_groups\[(\d+)\]\s*=\s*(.+?)\s*$')
    re_ref_marker = re.compile(r'^\s*refDB_marker\[(\d+)\]\s*=\s*(\d+)\s*$')
    re_ref_path = re.compile(r'^\s*refDB_path\[(\d+)\]\s*=\s*(.+?)\s*$')
    re_hie_marker = re.compile(r'^\s*hieDB_marker\[(\d+)\]\s*=\s*(\d+)\s*$')
    re_hie_path = re.compile(r'^\s*hieDB_path\[(\d+)\]\s*=\s*(.+?)\s*$')

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

            m = re_hie_marker.match(line)
            if m:
                hiedb_marker[int(m.group(1))] = int(m.group(2))
                continue

            m = re_hie_path.match(line)
            if m:
                hiedb_path[int(m.group(1))] = strip_quotes(m.group(2))
                continue

            m = re_old_marker.match(line)
            if m:
                old_marker = strip_quotes(m.group(1))
                continue

            m = re_old_ref.match(line)
            if m:
                old_refdbs[int(m.group(1))] = strip_quotes(m.group(2))
                continue

    def add_db_set(refdbs, marker_map, path_map, offset=0):
        all_indices = sorted(set(marker_map) | set(path_map))
        for src_i in all_indices:
            if src_i not in marker_map:
                raise SystemExit(f"ERROR: database marker[{src_i}] missing in config")
            if src_i not in path_map:
                raise SystemExit(f"ERROR: database path[{src_i}] missing in config")

            marker_i = marker_map[src_i]
            if marker_i not in marker_groups:
                raise SystemExit(
                    f"ERROR: database marker[{src_i}]={marker_i}, "
                    f"but marker_groups[{marker_i}] is not defined"
                )
            refdbs[marker_i][offset + src_i] = path_map[src_i]

    if marker_groups:
        refdbs = defaultdict(dict)

        if db_source in ("refDB", "both"):
            add_db_set(refdbs, refdb_marker, refdb_path, offset=0)

        if db_source in ("hieDB", "both"):
            # Offset avoids index collision if both refDB and hieDB are used.
            add_db_set(refdbs, hiedb_marker, hiedb_path, offset=100000)

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
        "ERROR: no marker_groups[] with selected database arrays found, "
        "and no old-style refDBs[] fallback found in config"
    )


def fasta_ids(fasta):
    ids = []

    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0].split(";")[0])

    return ids


def summarize_hits(hits):
    if not hits:
        return {
            "hits": 0,
            "unique_asvs": 0,
            "min_identity": None,
            "mean_identity": None,
            "max_identity": None,
        }

    unique_asvs = {h["asv_id"] for h in hits}
    identities = [h["identity"] for h in hits]

    return {
        "hits": len(hits),
        "unique_asvs": len(unique_asvs),
        "min_identity": min(identities),
        "mean_identity": sum(identities) / len(identities),
        "max_identity": max(identities),
    }


def run_vsearch(vsearch, query, db, out_uc, threshold, threads, log, strand="both", vsearch_extra=None):
    cmd = [
        vsearch,
        "--usearch_global", query,
        "--db", db,
        "--id", f"0.{threshold}",
        "--uc", out_uc,
        "--strand", strand,
        "--threads", str(threads),
    ]

    if vsearch_extra:
        cmd.extend(shlex.split(vsearch_extra))

    print("--- Starting marker assignment DB search", flush=True)
    print(f"    DB: {db}", flush=True)
    print(f"    UC output: {out_uc}", flush=True)
    print(f"    Log: {log}", flush=True)
    print("    Command: " + " ".join(shlex.quote(x) for x in cmd), flush=True)

    t0 = time.time()

    with open(log, "w") as err:
        err.write("Command: " + " ".join(shlex.quote(x) for x in cmd) + "\n")
        err.flush()

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

        if proc.stderr is None:
            raise RuntimeError("Failed to capture vsearch stderr")

        for line in proc.stderr:
            err.write(line)
            err.flush()

            lower = line.lower()
            # Print relevant progress/status lines live without flooding too much.
            if (
                "vsearch" in lower
                or "reading" in lower
                or "searching" in lower
                or "matching" in lower
                or "query" in lower
                or "hits" in lower
                or "warning" in lower
                or "error" in lower
                or "fatal" in lower
                or "%" in line
            ):
                print("    " + line.rstrip(), flush=True)

        return_code = proc.wait()

    elapsed = time.time() - t0

    if return_code != 0:
        print(f"--- DB search failed after {elapsed:.1f} sec", flush=True)
        raise subprocess.CalledProcessError(return_code, cmd)

    print(f"--- Finished DB search in {elapsed:.1f} sec", flush=True)
    return elapsed


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
    ap.add_argument(
        "--db-source",
        choices=["refDB", "hieDB", "both"],
        default="refDB",
        help="Database set to use for marker assignment. Use hieDB for one global DB per marker. Default: refDB."
    )
    ap.add_argument(
        "--strand",
        choices=["plus", "both"],
        default="both",
        help="VSEARCH search strand for marker assignment. 'plus' is faster if query and DB orientations are normalized. Default: both."
    )
    ap.add_argument(
        "--vsearch-extra",
        default="",
        help="Extra options passed to vsearch, e.g. '--qmask none --dbmask none --maxaccepts 1 --maxrejects 32'."
    )
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs("logs/marker_split", exist_ok=True)

    print("-- parsing marker assignment config", flush=True)
    markers, refdbs = parse_bash_config_arrays(args.config, db_source=args.db_source)

    if len(markers) < 1:
        raise SystemExit("ERROR: no marker groups found in config")

    for marker_i, marker_name in markers:
        if marker_i not in refdbs or not refdbs[marker_i]:
            raise SystemExit(
                f"ERROR: marker group {marker_i} ({marker_name}) has no database entries for db_source={args.db_source}"
            )

    print("-- loading ASV IDs", flush=True)
    all_asvs = fasta_ids(args.asv_fasta)
    all_hits = []

    print("-- marker assignment settings", flush=True)
    print(f"ASV FASTA: {args.asv_fasta}", flush=True)
    print(f"ASVs loaded: {len(all_asvs)}", flush=True)
    print(f"Config: {args.config}", flush=True)
    print(f"Threshold: {args.threshold}", flush=True)
    print(f"Min delta: {args.min_delta}", flush=True)
    print(f"Threads: {args.threads}", flush=True)
    print(f"DB source: {args.db_source}", flush=True)
    print(f"VSEARCH strand: {args.strand}", flush=True)
    if args.vsearch_extra:
        print(f"VSEARCH extra options: {args.vsearch_extra}", flush=True)

    print("-- marker groups detected", flush=True)
    for marker_i, marker_name in markers:
        print(f"marker_group[{marker_i}] = {marker_name}", flush=True)
        for db_i in sorted(refdbs.get(marker_i, {})):
            print(f"  DB {db_i}: {refdbs[marker_i][db_i]}", flush=True)

    for marker_i, marker_name in markers:
        marker_safe = safe_name(marker_name)
        marker_start_hits = len(all_hits)

        print(f"-- starting searches for marker group: {marker_name}", flush=True)

        for db_i in sorted(refdbs.get(marker_i, {})):
            db = refdbs[marker_i][db_i]

            uc = os.path.join(args.outdir, f"{marker_safe}.db{db_i}.uc")
            log = os.path.join("logs", "marker_split", f"{marker_safe}.db{db_i}.log")

            print(f"-- marker assignment search: marker={marker_name} db={db_i}", flush=True)
            print(f"DB: {db}", flush=True)

            run_vsearch(
                args.vsearch,
                args.asv_fasta,
                db,
                uc,
                args.threshold,
                args.threads,
                log,
                strand=args.strand,
                vsearch_extra=args.vsearch_extra,
            )

            print(f"--- Parsing UC file: {uc}", flush=True)
            db_hits = parse_uc_hits(uc, marker_name, db, db_i)
            all_hits.extend(db_hits)

            s = summarize_hits(db_hits)
            print(f"--- Parsed UC for marker={marker_name} db={db_i}", flush=True)
            print(f"    hit rows: {s['hits']}", flush=True)
            print(f"    unique ASVs with hit: {s['unique_asvs']} / {len(all_asvs)}", flush=True)
            if s["hits"] > 0:
                print(
                    "    identity min/mean/max: "
                    f"{s['min_identity']:.2f} / {s['mean_identity']:.2f} / {s['max_identity']:.2f}",
                    flush=True,
                )
            print(f"    cumulative hit rows so far: {len(all_hits)}", flush=True)

        marker_hits = all_hits[marker_start_hits:]
        marker_unique_asvs = {h["asv_id"] for h in marker_hits}
        print(
            f"-- marker {marker_name} cumulative unique ASVs with hits: "
            f"{len(marker_unique_asvs)} / {len(all_asvs)}",
            flush=True,
        )

    print("-- aggregating marker hits by ASV", flush=True)
    hits_by_asv = defaultdict(list)

    all_hits_path = os.path.join(args.outdir, "all_marker_hits.tsv")
    print(f"-- writing all marker hits: {all_hits_path}", flush=True)
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

    print("-- assigning best marker per ASV", flush=True)
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

        for asv_i, asv in enumerate(all_asvs, start=1):
            if asv_i % 10000 == 0:
                print(f"   assigned {asv_i} / {len(all_asvs)} ASVs", flush=True)

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

    print("-- writing marker-specific ID files", flush=True)
    for _, marker_name in markers:
        marker_safe = safe_name(marker_name)

        with open(os.path.join(args.outdir, f"{marker_safe}.ids"), "w") as out:
            for asv in assigned.get(marker_name, []):
                out.write(asv + "\n")

    for special in ["ambiguous", "unassigned"]:
        with open(os.path.join(args.outdir, f"{special}.ids"), "w") as out:
            for asv in assigned.get(special, []):
                out.write(asv + "\n")

    print("-- marker assignment summary", flush=True)
    total_assigned = 0

    for _, marker_name in markers:
        n = len(assigned.get(marker_name, []))
        total_assigned += n
        print(f"{marker_name}: {n}", flush=True)

    n_ambiguous = len(assigned.get("ambiguous", []))
    n_unassigned = len(assigned.get("unassigned", []))

    print(f"ambiguous: {n_ambiguous}", flush=True)
    print(f"unassigned: {n_unassigned}", flush=True)
    print(f"total ASVs: {len(all_asvs)}", flush=True)
    print(f"assigned to marker groups: {total_assigned}", flush=True)


if __name__ == "__main__":
    main()
