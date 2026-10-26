CRISPR SpCas9 Guide RNA Designer (CLI)
Scans FASTA or raw DNA sequences for 5'-NGG-3' PAM sites, extracts 20-nt
upstream protospacers, filters by GC%, and ranks candidates toward optimal GC%.
"""

import argparse
import csv
import sys
from pathlib import Path

# Default Benchmark: Human Beta-Globin (HBB) Coding Sequence (RefSeq NM_000518.4, 444 bp)
HBB_DEFAULT_SEQ = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC"
    "GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAG"
    "AGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAG"
    "GTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGAC"
    "AACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGAT"
    "CCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGC"
    "AAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAAT"
    "GCCCTGGCCCACAAGTATCACTAA"
)

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def parse_fasta(text: str) -> tuple[str, str]:
    """Parse FASTA string or raw sequence into (header, cleaned_dna_sequence)."""
    header = "NM_000518.4 Homo sapiens hemoglobin subunit beta (HBB), mRNA CDS"
    seq_lines = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip()
        else:
            seq_lines.append(line)

    raw_seq = "".join(seq_lines).upper().replace("U", "T")
    clean_seq = "".join(base for base in raw_seq if base in COMPLEMENT)
    return header, clean_seq


def reverse_complement(seq: str) -> str:
    """Return the 5'->3' reverse complement of a DNA string."""
    return "".join(COMPLEMENT.get(base, "N") for base in reversed(seq))


def calculate_gc(seq: str) -> float:
    """Return GC percentage of a DNA string."""
    if not seq:
        return 0.0
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100.0


def scan_strand(seq: str, strand: str, total_len: int) -> tuple[int, list[dict]]:
    """Scan a single DNA strand for 5'-NGG-3' PAM sites and 20-nt upstream guides."""
    pam_count = 0
    guides = []

    for i in range(len(seq) - 2):
        window = seq[i : i + 3]
        if window[1:] == "GG":
            pam_count += 1
            if i >= 20:
                guide_seq = seq[i - 20 : i]
                gc = calculate_gc(guide_seq)

                if strand == "+":
                    pam_pos_1based = i + 1
                    cut_site_1based = i - 2
                else:
                    pam_pos_1based = total_len - i
                    cut_site_1based = total_len - i + 3

                guides.append(
                    {
                        "strand": strand,
                        "scan_idx": i + 1,
                        "pam_pos": pam_pos_1based,
                        "cut_site": cut_site_1based,
                        "guide": guide_seq,
                        "pam": window,
                        "gc": gc,
                    }
                )

    return pam_count, guides


def run_pipeline(
    seq: str,
    strand_mode: str = "forward",
    min_gc: float = 40.0,
    max_gc: float = 60.0,
    opt_gc: float = 50.0,
    target_pos: int | None = 20,
):
    """Run full CRISPR discovery and ranking pipeline."""
    total_len = len(seq)
    fwd_pams, fwd_guides = scan_strand(seq, "+", total_len)
    rev_pams, rev_guides = scan_strand(reverse_complement(seq), "-", total_len)

    if strand_mode == "forward":
        total_pams = fwd_pams
        usable = fwd_guides
    elif strand_mode == "reverse":
        total_pams = rev_pams
        usable = rev_guides
    else:
        total_pams = fwd_pams + rev_pams
        usable = fwd_guides + rev_guides

    for g in usable:
        g["dist_from_opt"] = abs(g["gc"] - opt_gc)
        g["passes_gc"] = min_gc <= g["gc"] <= max_gc
        g["target_dist"] = (
            abs(g["cut_site"] - target_pos) if target_pos is not None else None
        )

    passed = [g for g in usable if g["passes_gc"]]
    ranked = sorted(
        passed,
        key=lambda x: (
            x["dist_from_opt"],
            x["scan_idx"],
            0 if x["strand"] == "+" else 1,
        ),
    )

    return total_pams, usable, ranked


def main():
    parser = argparse.ArgumentParser(
        description="CRISPR SpCas9 Guide RNA Designer (5'-NGG-3' PAM)"
    )
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        default=None,
        help="Path to input FASTA or TXT sequence file (defaults to HBB NM_000518.4)",
    )
    parser.add_argument(
        "-s",
        "--strand",
        choices=["forward", "reverse", "both"],
        default="forward",
        help="Strand scanning mode: forward (+), reverse (-), or both (default: forward)",
    )
    parser.add_argument(
        "--min-gc", type=float, default=40.0, help="Minimum GC%% (default: 40.0)"
    )
    parser.add_argument(
        "--max-gc", type=float, default=60.0, help="Maximum GC%% (default: 60.0)"
    )
    parser.add_argument(
        "--opt-gc", type=float, default=50.0, help="Optimal target GC%% (default: 50.0)"
    )
    parser.add_argument(
        "--target-pos",
        type=int,
        default=20,
        help="1-based target mutation coordinate to measure cut distance (default: 20)",
    )
    parser.add_argument(
        "--export",
        type=str,
        default=None,
        help="Optional output CSV filepath (e.g., guides.csv)",
    )

    args = parser.parse_args()

    if args.fasta:
        path = Path(args.fasta)
        if not path.is_file():
            print(f"Error: File not found: {args.fasta}", file=sys.stderr)
            sys.exit(1)
        header, seq = parse_fasta(path.read_text(encoding="utf-8"))
    else:
        header = "NM_000518.4 Homo sapiens hemoglobin subunit beta (HBB), mRNA CDS"
        seq = HBB_DEFAULT_SEQ

    total_pams, usable, ranked = run_pipeline(
        seq=seq,
        strand_mode=args.strand,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        opt_gc=args.opt_gc,
        target_pos=args.target_pos,
    )

    print("=" * 72)
    print(f"Target Sequence : {header}")
    print(f"Sequence Length : {len(seq)} bp ({calculate_gc(seq):.1f}% overall GC)")
    print(f"Strand Mode     : {args.strand.upper()}")
    print("=" * 72)
    print(f"Total PAM sites (NGG)                : {total_pams}")
    print(f"Usable guides (full 20bp upstream)   : {len(usable)}")
    print(
        f"Guides passing GC filter ({args.min_gc:.0f}%-{args.max_gc:.0f}%) : {len(ranked)}\n"
    )

    print(
        f"{'Rank':<5} {'Str':<4} {'PAM Pos':<8} {'Cut':<7} {'Dist':<7} {'20-nt Guide Sequence':<22} {'PAM':<5} {'GC %':<6}"
    )
    print("-" * 72)
    for idx, g in enumerate(ranked, 1):
        dist_str = f"{g['target_dist']}bp" if g["target_dist"] is not None else "N/A"
        print(
            f"{idx:<5} {g['strand']:<4} {g['pam_pos']:<8} {g['cut_site']:<7} {dist_str:<7} {g['guide']:<22} {g['pam']:<5} {g['gc']:5.1f}%"
        )

    if args.export:
        with open(args.export, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "Rank",
                    "Strand",
                    "PAM_Pos_1Based",
                    "Cut_Site_1Based",
                    "Dist_To_Target_Bp",
                    "Guide_Sequence_5to3",
                    "PAM",
                    "GC_Percent",
                ]
            )
            for idx, g in enumerate(ranked, 1):
                writer.writerow(
                    [
                        idx,
                        g["strand"],
                        g["pam_pos"],
                        g["cut_site"],
                        g["target_dist"] if g["target_dist"] is not None else "",
                        g["guide"],
                        g["pam"],
                        f"{g['gc']:.1f}",
                    ]
                )
        print(f"\n[+] Exported {len(ranked)} ranked guides to {args.export}")


if __name__ == "__main__":
    main()
