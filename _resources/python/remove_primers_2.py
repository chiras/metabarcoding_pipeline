#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import IUPACData
from tqdm import tqdm

# Primer presets. Values may be strings or lists of strings.
PRIMER_PRESETS = {
    "ITS2": {
        "fw": ["ATGCGATACTTGGTGTGAAT"],
        "rv": ["TCCTCCGCTTATTGATATGC"],
    },
    "fITS": { # FITS7 + RITS4
        "fw": ["GTGARTCATCGAATCTTTG"],
        "rv": ["TCCTCCGCTTATTGATATGC"],
    },
    "fITS+16S": {  # 16S: Kozich/EMP-style 515F + 806R; optional older 341F-like + 515R-like; ITS: fITS7 + ITS4/RITS4
        "fw": [
            "GTGARTCATCGAATCTTTG",    # fITS7
            "GTGCCAGCMGCCGCGGTA",     # 515F,  V4
            "CCTACGGGAGGCAGCAG",      # 341F-like, older/alternative 16S setup
        ],
        "rv": [
            "TCCTCCGCTTATTGATATGC",   # ITS4 / RITS4
            "GGACTACHVGGGTWTCTAAT",   # 806R
            "GTGCCAGCMGCCGCGGTAA",    # 515R-like, older/alternative 16S setup
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
    "16S": {  # accepted 16S primer variants: 341F/515R-like and 515F/806R V4
        "fw": [
            "CCTACGGGAGGCAGCAG",      # 341F-like
            "GTGCCAGCMGCCGCGGTA"      # 515F, as in Maroschek sheet
        ],
        "rv": [
            "GTGCCAGCMGCCGCGGTAA",    # 515R-like / reverse primer in older setup
            "GGACTACHVGGGTWTCTAAT"    # 806R
        ],
    }
}



def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove forward and reverse primers from FASTA sequences. "
                    "Supports IUPAC ambiguity codes, multiple primers, and optional terminal stagger bases."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument(
        "--fw",
        action="append",
        help="Forward primer sequence. Can be used multiple times or as comma-separated list.",
    )
    parser.add_argument(
        "--rv",
        action="append",
        help="Reverse primer sequence in original primer orientation. Can be used multiple times or as comma-separated list. "
             "The reverse complement is searched at the sequence end.",
    )
    parser.add_argument(
        "--marker",
        choices=PRIMER_PRESETS.keys(),
        help="Marker name to use preset forward and reverse primer sequences",
    )
    parser.add_argument(
        "--max-5p-prefix-stagger",
        type=int,
        default=0,
        help="Allow up to this many arbitrary bases before the forward primer and remove them together with the primer. "
             "Use this for heterogeneity spacers/stagger bases before FW primer.",
    )
    parser.add_argument(
        "--max-3p-suffix-stagger",
        type=int,
        default=0,
        help="Allow up to this many arbitrary bases after the reverse primer hit and remove them together with the primer. "
             "Use this for terminal bases after RV in merged reads.",
    )
    parser.add_argument(
        "--count-first",
        action="store_true",
        help="Count input records first to show a finite progress bar. Disabled by default to avoid a second full pass over large files.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1,
        help="Discard sequences shorter than this length after trimming. Default: 1.",
    )
    return parser.parse_args()


def flatten_primers(values: Optional[Iterable[str]]) -> List[str]:
    primers: List[str] = []
    if not values:
        return primers
    for value in values:
        for part in value.split(","):
            primer = part.strip().upper()
            if primer:
                primers.append(primer)
    # Deduplicate while preserving order
    seen = set()
    unique = []
    for primer in primers:
        if primer not in seen:
            seen.add(primer)
            unique.append(primer)
    return unique


def iupac_to_regex(iupac_seq: str) -> str:
    iupac_dict = dict(IUPACData.ambiguous_dna_values)
    iupac_dict["I"] = "ACGT"  # inosine in primers, treated as wildcard for trimming

    pieces = []
    for nuc in iupac_seq.upper():
        if nuc not in iupac_dict:
            raise ValueError(f"Unsupported base/IUPAC code in primer: {nuc!r}")
        pieces.append(f"[{iupac_dict[nuc]}]")
    return "".join(pieces)


def compile_fw_regex(fw_primers: List[str], max_prefix_stagger: int) -> Optional[re.Pattern]:
    if not fw_primers:
        return None
    alternatives = "|".join(iupac_to_regex(p) for p in fw_primers)
    # Remove optional heterogeneity spacer plus primer from the beginning.
    return re.compile(rf"^[ACGTN]{{0,{max_prefix_stagger}}}(?:{alternatives})", re.IGNORECASE)


def reverse_complement_primer(primer: str) -> str:
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "U": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
        "I": "I",  # inosine, later interpreted as [ACGT]
    }

    rc = []
    for nuc in reversed(primer.upper()):
        if nuc not in complement:
            raise ValueError(f"Unsupported base/IUPAC code in primer: {nuc!r}")
        rc.append(complement[nuc])
    return "".join(rc)


def compile_rv_regex(rv_primers: List[str], max_suffix_stagger: int) -> Optional[re.Pattern]:
    if not rv_primers:
        return None
    rv_rc_primers = [reverse_complement_primer(p) for p in rv_primers]
    alternatives = "|".join(iupac_to_regex(p) for p in rv_rc_primers)
    # Remove reverse-complement primer plus optional terminal heterogeneity spacer from the end.
    return re.compile(rf"(?:{alternatives})[ACGTN]{{0,{max_suffix_stagger}}}$", re.IGNORECASE)

def remove_primers(seq: Seq, fw_regex: Optional[re.Pattern], rv_regex: Optional[re.Pattern]) -> Tuple[Seq, dict]:
    trimmed = {"start": False, "end": False}
    seq_str = str(seq)

    if fw_regex:
        new_seq, n = fw_regex.subn("", seq_str, count=1)
        if n:
            seq_str = new_seq
            trimmed["start"] = True

    if rv_regex:
        new_seq, n = rv_regex.subn("", seq_str, count=1)
        if n:
            seq_str = new_seq
            trimmed["end"] = True

    return Seq(seq_str), trimmed


def get_primers(args) -> Tuple[List[str], List[str]]:
    if args.marker:
        preset = PRIMER_PRESETS[args.marker]
        fw_primers = flatten_primers(preset["fw"])
        rv_primers = flatten_primers(preset["rv"])
        # Allow user-provided primers to extend the preset.
        fw_primers.extend(p for p in flatten_primers(args.fw) if p not in fw_primers)
        rv_primers.extend(p for p in flatten_primers(args.rv) if p not in rv_primers)
    else:
        fw_primers = flatten_primers(args.fw)
        rv_primers = flatten_primers(args.rv)

    if not fw_primers and not rv_primers:
        raise ValueError("Provide at least one --fw/--rv primer or choose --marker.")
    return fw_primers, rv_primers


def main():
    args = parse_args()
    fw_primers, rv_primers = get_primers(args)

    fw_regex = compile_fw_regex(fw_primers, args.max_5p_prefix_stagger)
    rv_regex = compile_rv_regex(rv_primers, args.max_3p_suffix_stagger)

    total = None
    if args.count_first:
        total = sum(1 for _ in SeqIO.parse(args.input, "fasta"))

    trimmed_start_count = 0
    trimmed_end_count = 0
    trimmed_any_count = 0
    untrimmed_count = 0
    written_count = 0
    discarded_short_count = 0

    with open(args.output, "w") as out_handle:
        iterator = SeqIO.parse(args.input, "fasta")
        for record in tqdm(iterator, total=total, desc="Processing sequences", unit="seq"):
            new_seq, trimmed = remove_primers(record.seq, fw_regex, rv_regex)

            if len(new_seq) < args.min_length:
                discarded_short_count += 1
                continue

            if trimmed["start"]:
                trimmed_start_count += 1
            if trimmed["end"]:
                trimmed_end_count += 1
            if trimmed["start"] or trimmed["end"]:
                trimmed_any_count += 1
            else:
                untrimmed_count += 1

            new_record = SeqRecord(new_seq, id=record.id, name=record.name, description=record.description)
            SeqIO.write(new_record, out_handle, "fasta")
            written_count += 1

    print(f"Forward primers used: {', '.join(fw_primers) if fw_primers else 'none'}")
    print(f"Reverse primers used: {', '.join(rv_primers) if rv_primers else 'none'}")
    print(f"Sequences written: {written_count}")
    print(f"Sequences discarded because shorter than --min-length: {discarded_short_count}")
    print(f"Sequences trimmed at start: {trimmed_start_count}")
    print(f"Sequences trimmed at end: {trimmed_end_count}")
    print(f"Sequences trimmed at either end: {trimmed_any_count}")
    print(f"Sequences not trimmed: {untrimmed_count}")


if __name__ == "__main__":
    main()
