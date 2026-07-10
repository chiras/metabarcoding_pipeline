#!/usr/bin/env python3
"""Primer removal for merged metabarcoding FASTA files.

Key features
- IUPAC-aware primer matching
- optional heterogeneity/stagger bases before 5' primers and after 3' primers
- merged-read orientation handling
- reverse-complements reverse-orientation reads after trimming by default
- optional mismatch-tolerant primer matching, default 0 mismatches
- streaming multiprocessing for large FASTA files
"""

import argparse
import multiprocessing as mp
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm


PRIMER_PRESETS = {
    "ITS2": {
        "fw": ["ATGCGATACTTGGTGTGAAT"],
        "rv": ["TCCTCCGCTTATTGATATGC"],
    },
    "fITS": {
        "fw": ["GTGARTCATCGAATCTTTG"],
        "rv": ["TCCTCCGCTTATTGATATGC"],
    },
    "fITS+16S": {
        "fw": [
            "GTGARTCATCGAATCTTTG",    # fITS7
            "GTGCCAGCMGCCGCGGTA",     # 515F, V4
            "CCTACGGGAGGCAGCAG",      # 341F-like, older/alternative 16S setup
        ],
        "rv": [
            "TCCTCCGCTTATTGATATGC",   # ITS4 / RITS4
            "GGACTACHVGGGTWTCTAAT",   # 806R
            "GTGCCAGCMGCCGCGGTAA",    # 515R-like, older/alternative 16S setup
        ],
    },
    "16S": {
        "fw": [
            "CCTACGGGAGGCAGCAG",      # 341F-like
            "GTGCCAGCMGCCGCGGTA",     # 515F, V4
        ],
        "rv": [
            "GTGCCAGCMGCCGCGGTAA",    # 515R-like / reverse primer in older setup
            "GGACTACHVGGGTWTCTAAT",   # 806R
        ],
    },
    "COI": {
        "fw": ["GGWACWGGWTGAACWGTWTAYCCYCC"],
        "rv": [
            "TAIACYTCIGGRTGICCRAARAAYCA",
            "TAAACTTCAGGGTGACCAAARAAYCA",
        ],
    },
    "COI-5P": {
        "fw": ["GGWACWGGWTGAACWGTWTAYCCYCC"],
        "rv": [
            "TAIACYTCIGGRTGICCRAARAAYCA",
            "TAAACTTCAGGGTGACCAAARAAYCA",
        ],
    },
}

IUPAC: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
    "I": "ACGT",  # inosine in primers, treated as wildcard
}

RC_TABLE = str.maketrans(
    "ACGTURYKMSWBDHVNIacgturykmswbdhvni",
    "TGCAAYRMKSWVHDBNItgcaayrmkswvhdbni",
)

WORKER_CONTEXT = None


@dataclass
class PrimerHit:
    start: int
    end: int
    mismatches: int
    primer_length: int


@dataclass
class TrimResult:
    seq: str
    orientation: str
    start_trimmed: bool
    end_trimmed: bool
    mismatches: int

    @property
    def score(self) -> int:
        return int(self.start_trimmed) + int(self.end_trimmed)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Remove primer sequences from merged FASTA reads. Supports IUPAC ambiguity codes, "
            "multiple primer pairs, heterogeneity spacers/stagger bases, both merged-read orientations, "
            "optional orientation normalization by reverse-complementing, and optional mismatch tolerance."
        )
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--fw", action="append", help="Forward primer. Can be repeated or comma-separated.")
    parser.add_argument(
        "--rv",
        action="append",
        help=(
            "Reverse primer in original primer orientation. Can be repeated or comma-separated. "
            "For forward-orientation merged reads, its reverse complement is searched at the sequence end."
        ),
    )
    parser.add_argument("--marker", choices=PRIMER_PRESETS.keys(), help="Primer preset marker name")
    parser.add_argument(
        "--max-5p-prefix-stagger",
        type=int,
        default=0,
        help="Allow and remove up to this many arbitrary bases before a 5' primer.",
    )
    parser.add_argument(
        "--max-3p-suffix-stagger",
        type=int,
        default=0,
        help="Allow and remove up to this many arbitrary bases after a 3' primer.",
    )
    parser.add_argument(
        "--max-primer-mismatches",
        type=int,
        default=0,
        help="Allowed mismatches per primer. Default 0. With 0, a fast exact IUPAC regex path is used.",
    )
    parser.add_argument(
        "--no-orient-normalize",
        action="store_true",
        help="Do not reverse-complement reverse-orientation reads after trimming. Default is to normalize orientation.",
    )
    parser.add_argument(
        "--discard-untrimmed",
        action="store_true",
        help="Discard reads where no primer evidence is found in either orientation.",
    )
    parser.add_argument("--min-length", type=int, default=1, help="Discard reads shorter than this after trimming.")
    parser.add_argument("--count-first", action="store_true", help="Count records first for finite progress bar.")
    parser.add_argument("--threads", type=int, default=1, help="Number of worker processes. Default 1.")
    parser.add_argument("--batch-size", type=int, default=5000, help="Records per multiprocessing batch. Default 5000.")
    return parser.parse_args()


def flatten_primers(values: Optional[Iterable[str]]) -> List[str]:
    primers: List[str] = []
    if not values:
        return primers
    for value in values:
        for part in value.split(","):
            primer = part.strip().upper().replace("U", "T")
            if primer:
                primers.append(primer)
    seen = set()
    unique = []
    for primer in primers:
        if primer not in seen:
            validate_primer(primer)
            seen.add(primer)
            unique.append(primer)
    return unique


def validate_primer(primer: str) -> None:
    for nuc in primer.upper():
        if nuc not in IUPAC:
            raise ValueError(f"Unsupported base/IUPAC code in primer {primer!r}: {nuc!r}")


def reverse_complement_primer(primer: str) -> str:
    return primer.translate(RC_TABLE)[::-1].upper().replace("U", "T")


def reverse_complement_sequence(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper().replace("U", "T")


def iupac_to_regex(iupac_seq: str) -> str:
    return "".join(f"[{IUPAC[nuc]}]" for nuc in iupac_seq.upper())


def compile_5p_regex(primers: Sequence[str], max_prefix_stagger: int) -> Optional[re.Pattern]:
    if not primers:
        return None
    alternatives = "|".join(iupac_to_regex(p) for p in primers)
    return re.compile(rf"^[ACGTN]{{0,{max_prefix_stagger}}}(?:{alternatives})", re.IGNORECASE)


def compile_3p_regex(primers: Sequence[str], max_suffix_stagger: int) -> Optional[re.Pattern]:
    if not primers:
        return None
    alternatives = "|".join(iupac_to_regex(p) for p in primers)
    return re.compile(rf"(?:{alternatives})[ACGTN]{{0,{max_suffix_stagger}}}$", re.IGNORECASE)


def regex_hit_5p(seq: str, pattern: Optional[re.Pattern]) -> Optional[PrimerHit]:
    if pattern is None:
        return None
    m = pattern.search(seq)
    if not m:
        return None
    return PrimerHit(start=m.start(), end=m.end(), mismatches=0, primer_length=m.end() - m.start())


def regex_hit_3p(seq: str, pattern: Optional[re.Pattern]) -> Optional[PrimerHit]:
    if pattern is None:
        return None
    m = pattern.search(seq)
    if not m:
        return None
    return PrimerHit(start=m.start(), end=m.end(), mismatches=0, primer_length=m.end() - m.start())


def trim_orientation_regex(seq: str, start_re: Optional[re.Pattern], end_re: Optional[re.Pattern], orientation: str) -> TrimResult:
    seq_str = seq.upper().replace("U", "T")
    start_hit = regex_hit_5p(seq_str, start_re)
    end_hit = regex_hit_3p(seq_str, end_re)

    left = start_hit.end if start_hit else 0
    right = end_hit.start if end_hit else len(seq_str)

    if right < left:
        # Conflicting hits; keep only the better-supported start hit to avoid negative sequence.
        right = len(seq_str)
        end_hit = None

    return TrimResult(
        seq=seq_str[left:right],
        orientation=orientation,
        start_trimmed=start_hit is not None,
        end_trimmed=end_hit is not None,
        mismatches=0,
    )


def iupac_base_match(seq_base: str, primer_code: str) -> bool:
    seq_base = seq_base.upper().replace("U", "T")
    primer_code = primer_code.upper().replace("U", "T")
    if seq_base == "N":
        return True
    return seq_base in IUPAC[primer_code]


def primer_match(seq_segment: str, primer: str, max_mismatches: int) -> Optional[int]:
    if len(seq_segment) != len(primer):
        return None
    mismatches = 0
    for sb, pc in zip(seq_segment.upper(), primer.upper()):
        if not iupac_base_match(sb, pc):
            mismatches += 1
            if mismatches > max_mismatches:
                return None
    return mismatches


def better_hit(candidate: PrimerHit, current: Optional[PrimerHit]) -> PrimerHit:
    if current is None:
        return candidate
    # Prefer fewer mismatches, then longer primer, then less stagger/closer to terminal edge.
    cand_key = (candidate.mismatches, -candidate.primer_length, candidate.start)
    curr_key = (current.mismatches, -current.primer_length, current.start)
    return candidate if cand_key < curr_key else current


def find_5p_primer(seq: str, primers: Sequence[str], max_prefix_stagger: int, max_mismatches: int) -> Optional[PrimerHit]:
    best: Optional[PrimerHit] = None
    n = len(seq)
    for stagger in range(max_prefix_stagger + 1):
        for primer in primers:
            plen = len(primer)
            if stagger + plen > n:
                continue
            segment = seq[stagger:stagger + plen]
            mm = primer_match(segment, primer, max_mismatches)
            if mm is not None:
                best = better_hit(PrimerHit(0, stagger + plen, mm, plen), best)
    return best


def find_3p_primer(seq: str, primers: Sequence[str], max_suffix_stagger: int, max_mismatches: int) -> Optional[PrimerHit]:
    best: Optional[PrimerHit] = None
    n = len(seq)
    for suffix in range(max_suffix_stagger + 1):
        for primer in primers:
            plen = len(primer)
            start = n - suffix - plen
            end = n - suffix
            if start < 0:
                continue
            segment = seq[start:end]
            mm = primer_match(segment, primer, max_mismatches)
            if mm is not None:
                # For 3' hits, prefer fewer mismatches/longer primer/less terminal suffix.
                candidate = PrimerHit(start, n, mm, plen)
                if best is None:
                    best = candidate
                else:
                    cand_key = (candidate.mismatches, -candidate.primer_length, n - candidate.end)
                    best_key = (best.mismatches, -best.primer_length, n - best.end)
                    if cand_key < best_key:
                        best = candidate
    return best


def trim_orientation_slow(
    seq: str,
    start_primers: Sequence[str],
    end_primers: Sequence[str],
    orientation: str,
    max_prefix_stagger: int,
    max_suffix_stagger: int,
    max_mismatches: int,
) -> TrimResult:
    seq_str = seq.upper().replace("U", "T")
    start_hit = find_5p_primer(seq_str, start_primers, max_prefix_stagger, max_mismatches)
    end_hit = find_3p_primer(seq_str, end_primers, max_suffix_stagger, max_mismatches)

    left = start_hit.end if start_hit else 0
    right = end_hit.start if end_hit else len(seq_str)

    if right < left:
        right = len(seq_str)
        end_hit = None

    return TrimResult(
        seq=seq_str[left:right],
        orientation=orientation,
        start_trimmed=start_hit is not None,
        end_trimmed=end_hit is not None,
        mismatches=(start_hit.mismatches if start_hit else 0) + (end_hit.mismatches if end_hit else 0),
    )


def choose_best_result(forward: TrimResult, reverse: TrimResult) -> Tuple[TrimResult, bool]:
    """Return best trim result and whether evidence was ambiguous."""
    if forward.score > reverse.score:
        return forward, False
    if reverse.score > forward.score:
        return reverse, False

    if forward.score == 0 and reverse.score == 0:
        return forward, False

    # Same number of trimmed ends; prefer fewer mismatches.
    if forward.mismatches < reverse.mismatches:
        return forward, False
    if reverse.mismatches < forward.mismatches:
        return reverse, False

    # Same evidence. Prefer the result with longer retained insert.
    if len(forward.seq) > len(reverse.seq):
        return forward, True
    if len(reverse.seq) > len(forward.seq):
        return reverse, True

    # Stable tie-breaker: forward.
    return forward, True


def remove_primers_fast_regex(seq: str, context: dict) -> Tuple[str, dict]:
    forward = trim_orientation_regex(
        seq,
        context["fw_start_re"],
        context["rv_rc_end_re"],
        "forward",
    )
    reverse = trim_orientation_regex(
        seq,
        context["rv_start_re"],
        context["fw_rc_end_re"],
        "reverse",
    )
    best, ambiguous = choose_best_result(forward, reverse)

    out_seq = best.seq
    reverse_complemented = False
    if best.orientation == "reverse" and context["orient_normalize"]:
        out_seq = reverse_complement_sequence(out_seq)
        reverse_complemented = True

    return out_seq, {
        "start": best.start_trimmed,
        "end": best.end_trimmed,
        "any": best.start_trimmed or best.end_trimmed,
        "orientation": best.orientation,
        "reverse_complemented": reverse_complemented,
        "ambiguous": ambiguous,
        "mismatches": best.mismatches,
    }


def remove_primers_slow_mismatch(seq: str, context: dict) -> Tuple[str, dict]:
    fw_primers = context["fw_primers"]
    rv_primers = context["rv_primers"]
    fw_rc_primers = context["fw_rc_primers"]
    rv_rc_primers = context["rv_rc_primers"]

    forward = trim_orientation_slow(
        seq,
        fw_primers,
        rv_rc_primers,
        "forward",
        context["max_5p_prefix_stagger"],
        context["max_3p_suffix_stagger"],
        context["max_primer_mismatches"],
    )
    reverse = trim_orientation_slow(
        seq,
        rv_primers,
        fw_rc_primers,
        "reverse",
        context["max_5p_prefix_stagger"],
        context["max_3p_suffix_stagger"],
        context["max_primer_mismatches"],
    )
    best, ambiguous = choose_best_result(forward, reverse)

    out_seq = best.seq
    reverse_complemented = False
    if best.orientation == "reverse" and context["orient_normalize"]:
        out_seq = reverse_complement_sequence(out_seq)
        reverse_complemented = True

    return out_seq, {
        "start": best.start_trimmed,
        "end": best.end_trimmed,
        "any": best.start_trimmed or best.end_trimmed,
        "orientation": best.orientation,
        "reverse_complemented": reverse_complemented,
        "ambiguous": ambiguous,
        "mismatches": best.mismatches,
    }


def empty_stats() -> dict:
    return {
        "trimmed_start_count": 0,
        "trimmed_end_count": 0,
        "trimmed_any_count": 0,
        "untrimmed_count": 0,
        "written_count": 0,
        "discarded_short_count": 0,
        "discarded_untrimmed_count": 0,
        "forward_orientation_count": 0,
        "reverse_orientation_count": 0,
        "reverse_complemented_count": 0,
        "ambiguous_orientation_count": 0,
        "total_primer_mismatches": 0,
    }


def add_stats(dest: dict, src: dict) -> None:
    for key, value in src.items():
        dest[key] = dest.get(key, 0) + value


def fasta_text(title: str, seq: str) -> str:
    return f">{title}\n{seq}\n"


def process_record(title: str, seq: str, context: dict) -> Tuple[Optional[str], dict]:
    stats = empty_stats()

    if context["use_fast_regex"]:
        new_seq, trimmed = remove_primers_fast_regex(seq, context)
    else:
        new_seq, trimmed = remove_primers_slow_mismatch(seq, context)

    if trimmed["orientation"] == "forward":
        stats["forward_orientation_count"] += 1
    elif trimmed["orientation"] == "reverse":
        stats["reverse_orientation_count"] += 1

    if trimmed["reverse_complemented"]:
        stats["reverse_complemented_count"] += 1
    if trimmed["ambiguous"]:
        stats["ambiguous_orientation_count"] += 1
    if trimmed["start"]:
        stats["trimmed_start_count"] += 1
    if trimmed["end"]:
        stats["trimmed_end_count"] += 1
    if trimmed["any"]:
        stats["trimmed_any_count"] += 1
    else:
        stats["untrimmed_count"] += 1
    stats["total_primer_mismatches"] += trimmed["mismatches"]

    if context["discard_untrimmed"] and not trimmed["any"]:
        stats["discarded_untrimmed_count"] += 1
        return None, stats

    if len(new_seq) < context["min_length"]:
        stats["discarded_short_count"] += 1
        return None, stats

    stats["written_count"] += 1
    return fasta_text(title, new_seq), stats


def init_worker(context: dict) -> None:
    global WORKER_CONTEXT
    ctx = dict(context)

    fw_primers = ctx["fw_primers"]
    rv_primers = ctx["rv_primers"]
    fw_rc_primers = [reverse_complement_primer(p) for p in fw_primers]
    rv_rc_primers = [reverse_complement_primer(p) for p in rv_primers]

    ctx["fw_rc_primers"] = fw_rc_primers
    ctx["rv_rc_primers"] = rv_rc_primers

    if ctx["use_fast_regex"]:
        ctx["fw_start_re"] = compile_5p_regex(fw_primers, ctx["max_5p_prefix_stagger"])
        ctx["rv_start_re"] = compile_5p_regex(rv_primers, ctx["max_5p_prefix_stagger"])
        ctx["rv_rc_end_re"] = compile_3p_regex(rv_rc_primers, ctx["max_3p_suffix_stagger"])
        ctx["fw_rc_end_re"] = compile_3p_regex(fw_rc_primers, ctx["max_3p_suffix_stagger"])

    WORKER_CONTEXT = ctx


def process_batch(batch: List[Tuple[str, str]]) -> Tuple[str, dict, int]:
    out_parts: List[str] = []
    batch_stats = empty_stats()
    for title, seq in batch:
        out_text, rec_stats = process_record(title, seq, WORKER_CONTEXT)
        add_stats(batch_stats, rec_stats)
        if out_text is not None:
            out_parts.append(out_text)
    return "".join(out_parts), batch_stats, len(batch)


def batched_fasta_records(handle, batch_size: int):
    batch: List[Tuple[str, str]] = []
    for title, seq in SimpleFastaParser(handle):
        batch.append((title, seq))
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def get_primers(args) -> Tuple[List[str], List[str]]:
    if args.marker:
        preset = PRIMER_PRESETS[args.marker]
        fw_primers = flatten_primers(preset["fw"])
        rv_primers = flatten_primers(preset["rv"])
        for primer in flatten_primers(args.fw):
            if primer not in fw_primers:
                fw_primers.append(primer)
        for primer in flatten_primers(args.rv):
            if primer not in rv_primers:
                rv_primers.append(primer)
    else:
        fw_primers = flatten_primers(args.fw)
        rv_primers = flatten_primers(args.rv)

    if not fw_primers and not rv_primers:
        raise ValueError("Provide at least one --fw/--rv primer or choose --marker.")
    return fw_primers, rv_primers


def main() -> None:
    args = parse_args()

    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.batch_size < 1:
        raise ValueError("--batch-size must be >= 1")
    if args.max_5p_prefix_stagger < 0 or args.max_3p_suffix_stagger < 0:
        raise ValueError("stagger values must be >= 0")
    if args.max_primer_mismatches < 0:
        raise ValueError("--max-primer-mismatches must be >= 0")

    fw_primers, rv_primers = get_primers(args)
    use_fast_regex = args.max_primer_mismatches == 0
    orient_normalize = not args.no_orient_normalize

    total = None
    if args.count_first:
        with open(args.input) as in_handle:
            total = sum(1 for _title, _seq in SimpleFastaParser(in_handle))

    worker_context = {
        "fw_primers": fw_primers,
        "rv_primers": rv_primers,
        "max_5p_prefix_stagger": args.max_5p_prefix_stagger,
        "max_3p_suffix_stagger": args.max_3p_suffix_stagger,
        "max_primer_mismatches": args.max_primer_mismatches,
        "orient_normalize": orient_normalize,
        "use_fast_regex": use_fast_regex,
        "discard_untrimmed": args.discard_untrimmed,
        "min_length": args.min_length,
    }

    stats = empty_stats()

    with open(args.input) as in_handle, open(args.output, "w") as out_handle:
        if args.threads == 1:
            init_worker(worker_context)
            iterator = SimpleFastaParser(in_handle)
            progress = tqdm(iterator, total=total, desc="Processing sequences", unit="seq")
            for title, seq in progress:
                out_text, rec_stats = process_record(title, seq, WORKER_CONTEXT)
                add_stats(stats, rec_stats)
                if out_text is not None:
                    out_handle.write(out_text)
        else:
            batches = batched_fasta_records(in_handle, args.batch_size)
            with mp.Pool(processes=args.threads, initializer=init_worker, initargs=(worker_context,)) as pool:
                progress = tqdm(total=total, desc="Processing sequences", unit="seq")
                for out_text, batch_stats, batch_n in pool.imap(process_batch, batches, chunksize=1):
                    out_handle.write(out_text)
                    add_stats(stats, batch_stats)
                    progress.update(batch_n)
                progress.close()

    print(f"Forward primers used: {', '.join(fw_primers) if fw_primers else 'none'}")
    print(f"Reverse primers used: {', '.join(rv_primers) if rv_primers else 'none'}")
    print(f"Max primer mismatches allowed per primer: {args.max_primer_mismatches}")
    print(f"Matching mode: {'fast exact IUPAC regex' if use_fast_regex else 'slow mismatch-tolerant scan'}")
    print(f"Worker processes: {args.threads}")
    print(f"Batch size: {args.batch_size}")
    print(f"Orientation normalization by reverse-complementing reverse-orientation reads: {'yes' if orient_normalize else 'no'}")
    print(f"Sequences written: {stats['written_count']}")
    print(f"Sequences discarded because shorter than --min-length: {stats['discarded_short_count']}")
    print(f"Sequences discarded because untrimmed: {stats['discarded_untrimmed_count']}")
    print(f"Sequences assigned forward orientation: {stats['forward_orientation_count']}")
    print(f"Sequences assigned reverse orientation: {stats['reverse_orientation_count']}")
    print(f"Sequences reverse-complemented: {stats['reverse_complemented_count']}")
    print(f"Sequences with ambiguous orientation evidence: {stats['ambiguous_orientation_count']}")
    print(f"Sequences trimmed at start: {stats['trimmed_start_count']}")
    print(f"Sequences trimmed at end: {stats['trimmed_end_count']}")
    print(f"Sequences trimmed at either end: {stats['trimmed_any_count']}")
    print(f"Sequences not trimmed: {stats['untrimmed_count']}")
    print(f"Total primer mismatches accepted: {stats['total_primer_mismatches']}")


if __name__ == "__main__":
    main()
