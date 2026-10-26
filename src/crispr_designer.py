CRISPR SpCas9 Guide RNA Designer (Unified Engine)
=================================================
100% parity with the interactive web application.
Scans forward (+) and reverse (-) strands for SpCas9 PAM sites (NGG / CCN),
extracts 20-nt spacers, maps 1-based forward-strand coordinates and Cas9 cut
sites (-3 bp from PAM), calculates distance to a target locus, and filters/ranks
candidates by GC content and target proximity.
"""

import argparse
import csv
import sys
from pathlib import Path

# Default Benchmark: Human Beta-Globin (HBB) Coding Sequence (RefSeq NM_000518.4)
DEFAULT_HBB_SEQ = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC"
    "GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAG"
    "AGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAG"
    "GTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGAC"
    "AACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGAT"
    "CCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGC"
    "AAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAAT"
    "GCCCTGGCCCACAAGTATCACTAA"
)


def clean_sequence(raw_text: str) -> str:
    """Strip FASTA headers and non-ATCG characters."""
    lines = [
        line.strip()
        for line in raw_text.strip().splitlines()
        if not line.strip().startswith(">")
    ]
    seq = "".join(lines).upper()
    return "".join(b for b in seq if b in "ATCG")


def reverse_complement(seq: str) -> str:
    """Return the 5' -> 3' reverse complement of a DNA string."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))


def calculate_gc(seq: str) -> float:
    """Return GC percentage of a DNA string."""
    if not seq:
        return 0.0
    return ((seq.count("G") + seq.count("C")) / len(seq)) * 100.0


def scan_crispr_guides(
    raw_seq: str,
    strand_mode: str = "both",
    min_gc: float = 40.0,
    max_gc: float = 60.0,
    ideal_gc: float = 50.0,
    target_pos: int | None = 20,
    sort_by: str = "gc",
) -> dict:
    """
    Scan a DNA sequence for SpCas9 guides with identical coordinate math to the Web App.

    Coordinate Conventions (1-based on forward strand):
    - Forward (+) NGG at 0-based index `i`:
        PAM = seq[i:i+3] (1-based: i+1 .. i+3)
        Guide (5'->3') = seq[i-20:i] (1-based: i-19 .. i)
        Cut site = i - 2 (3 bp upstream of PAM, i.e., after base i-2)
    - Reverse (-) CCN on forward strand at 0-based index `i`:
        Forward CCN = seq[i:i+3] -> Reverse PAM (5'->3') = rev_comp(seq[i:i+3])
        Forward protospacer = seq[i+3:i+23] (1-based: i+4 .. i+23)
        Guide (5'->3') = rev_comp(seq[i+3:i+23])
        Cut site = i + 6 (3 bp upstream of PAM on reverse strand, mapped to forward 1-based base)
    """
    seq = clean_sequence(raw_seq)
    seq_len = len(seq)

    stats = {
        "seq_length": seq_len,
        "fwd_pams": 0,
        "rev_pams": 0,
        "fwd_usable": 0,
        "rev_usable": 0,
        "fwd_passed": 0,
        "rev_passed": 0,
    }
    usable_guides = []

    # 1. Forward Strand (+) Scan (NGG)
    if strand_mode in ("both", "+", "forward"):
        for i in range(seq_len - 2):
            window = seq[i : i + 3]
            if window[1:] == "GG":
                stats["fwd_pams"] += 1
                if i >= 20:
                    stats["fwd_usable"] += 1
                    guide_seq = seq[i - 20 : i]
                    gc = calculate_gc(guide_seq)
                    pam_pos_1b = i + 1
                    guide_start_1b = i - 19
                    guide_end_1b = i
                    cut_site_1b = i - 2  # 3 bp upstream of PAM start (pam_pos_1b - 3)
                    dist_target = (
                        abs(cut_site_1b - target_pos) if target_pos is not None else None
                    )
                    passes_gc = min_gc <= gc <= max_gc
                    if passes_gc:
                        stats["fwd_passed"] += 1

                    usable_guides.append(
                        {
                            "strand": "+",
                            "pam_pos": pam_pos_1b,
                            "guide_start": guide_start_1b,
                            "guide_end": guide_end_1b,
                            "cut_site": cut_site_1b,
                            "guide": guide_seq,
                            "pam": window,
                            "gc": round(gc, 2),
                            "dist_from_ideal_gc": round(abs(gc - ideal_gc), 2),
                            "dist_from_target": dist_target,
                            "passes_gc": passes_gc,
                        }
                    )

    # 2. Reverse Strand (-) Scan (CCN on forward strand)
    if strand_mode in ("both", "-", "reverse"):
        for i in range(seq_len - 2):
            window = seq[i : i + 3]
            if window[:2] == "CC":
                stats["rev_pams"] += 1
                if i + 23 <= seq_len:
                    stats["rev_usable"] += 1
                    fwd_protospacer = seq[i + 3 : i + 23]
                    guide_seq = reverse_complement(fwd_protospacer)
                    rev_pam = reverse_complement(window)
                    gc = calculate_gc(guide_seq)
                    pam_pos_1b = i + 1  # 1-based start of CCN motif on forward strand
                    guide_start_1b = i + 4
                    guide_end_1b = i + 23
                    cut_site_1b = i + 6  # 3 bp inside protospacer from PAM
                    dist_target = (
                        abs(cut_site_1b - target_pos) if target_pos is not None else None
                    )
                    passes_gc = min_gc <= gc <= max_gc
                    if passes_gc:
                        stats["rev_passed"] += 1

                    usable_guides.append(
                        {
                            "strand": "-",
                            "pam_pos": pam_pos_1b,
                            "guide_start": guide_start_1b,
                            "guide_end": guide_end_1b,
                            "cut_site": cut_site_1b,
                            "guide": guide_seq,
                            "pam": rev_pam,
                            "gc": round(gc, 2),
                            "dist_from_ideal_gc": round(abs(gc - ideal_gc), 2),
                            "dist_from_target": dist_target,
                            "passes_gc": passes_gc,
                        }
                    )

    # Filter by GC window
    passed_guides = [g for g in usable_guides if g["passes_gc"]]

    # Sort candidates identically to Web App
    if sort_by == "target" and target_pos is not None:
        passed_guides.sort(
            key=lambda x: (
                x["dist_from_target"],
                x["dist_from_ideal_gc"],
                x["pam_pos"],
            )
        )
    elif sort_by == "position":
        passed_guides.sort(key=lambda x: (x["pam_pos"], x["dist_from_ideal_gc"]))
    else:  # Default: sort by proximity to ideal 50% GC, then position
        passed_guides.sort(
            key=lambda x: (
                x["dist_from_ideal_gc"],
                x["pam_pos"],
                0 if x["strand"] == "+" else 1,
            )
        )

    stats["total_pams"] = stats["fwd_pams"] + stats["rev_pams"]
    stats["total_usable"] = stats["fwd_usable"] + stats["rev_usable"]
    stats["total_passed"] = stats["fwd_passed"] + stats["rev_passed"]

    return {
        "stats": stats,
        "usable_guides": usable_guides,
        "passed_guides": passed_guides,
    }


def print_report(result: dict, strand_mode: str, min_gc: float, max_gc: float, target_pos: int | None):
    s = result["stats"]
    passed = result["passed_guides"]

    print("=" * 84)
    print(f"CRISPR SpCas9 GUIDE RNA REPORT (Strand Mode: {strand_mode.upper()}, GC: {min_gc:.0f}%-{max_gc:.0f}%)")
    print("=" * 84)
    print(f"Sequence Length                  : {s['seq_length']} bp")
    print(f"Forward (+) PAMs / Usable / Pass : {s['fwd_pams']} PAMs | {s['fwd_usable']} Usable | {s['fwd_passed']} Passed")
    print(f"Reverse (-) PAMs / Usable / Pass : {s['rev_pams']} PAMs | {s['rev_usable']} Usable | {s['rev_passed']} Passed")
    print(f"Active Selection Totals          : {s['total_pams']} PAMs | {s['total_usable']} Usable | {s['total_passed']} Passed")
    print("-" * 84)

    header = (
        f"{'Rank':<5} {'Str':<4} {'Span (1b)':<12} {'PAM Pos':<8} "
        f"{'Cut':<6} {'ΔTarget':<8} {'20-nt Guide (5\'->3\')':<22} {'PAM':<5} {'GC%':<6}"
    )
    print(header)
    print("-" * 84)

    for rank, g in enumerate(passed, 1):
        span_str = f"{g['guide_start']}-{g['guide_end']}"
        dist_str = f"{g['dist_from_target']} bp" if g["dist_from_target"] is not None else "N/A"
        print(
            f"{rank:<5} {g['strand']:<4} {span_str:<12} {g['pam_pos']:<8} "
            f"{g['cut_site']:<6} {dist_str:<8} {g['guide']:<22} {g['pam']:<5} {g['gc']:5.1f}%"
        )


def main():
    parser = argparse.ArgumentParser(description="CRISPR SpCas9 Guide RNA Designer (Unified Engine)")
    parser.add_argument("-f", "--fasta", type=str, help="Path to input FASTA or text file (defaults to HBB NM_000518.4)")
    parser.add_argument("-s", "--strand", choices=["both", "+", "-"], default="both", help="Strand to scan: 'both', '+', or '-' (default: both)")
    parser.add_argument("--min-gc", type=float, default=40.0, help="Minimum GC percentage (default: 40.0)")
    parser.add_argument("--max-gc", type=float, default=60.0, help="Maximum GC percentage (default: 60.0)")
    parser.add_argument("-t", "--target", type=int, default=20, help="Target base position (1-based) for distance calculation (default: 20 for HBB E6V)")
    parser.add_argument("--sort", choices=["gc", "target", "position"], default="gc", help="Ranking criterion: 'gc' (closest to 50%%), 'target' (closest cut site), or 'position'")
    parser.add_argument("-o", "--out-csv", type=str, help="Optional path to export ranked passing guides as CSV")

    args = parser.parse_args()

    if args.fasta:
        raw_seq = Path(args.fasta).read_text(encoding="utf-8")
    else:
        raw_seq = DEFAULT_HBB_SEQ

    result = scan_crispr_guides(
        raw_seq=raw_seq,
        strand_mode=args.strand,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        target_pos=args.target,
        sort_by=args.sort,
    )

    print_report(result, args.strand, args.min_gc, args.max_gc, args.target)

    if args.out_csv:
        with open(args.out_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "strand",
                    "guide_start",
                    "guide_end",
                    "pam_pos",
                    "cut_site",
                    "dist_from_target",
                    "guide",
                    "pam",
                    "gc",
                    "dist_from_ideal_gc",
                ],
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(result["passed_guides"])
        print(f"\n[+] Exported {len(result['passed_guides'])} guides to {args.out_csv}")


if __name__ == "__main__":
    main()
