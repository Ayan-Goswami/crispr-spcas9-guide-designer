# Project 03: The CRISPR Guide Designer
# Target: Human HBB coding sequence (RefSeq NM_000518.4)

hbb_sequence = (
    "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGA"
    "GGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGC"
    "AGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATG"
    "CTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGC"
    "TCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGAT"
    "CCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCA"
    "CCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCA"
    "CTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACT"
    "GGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGC"
)


def calculate_gc(seq):
    """Return GC percentage of a DNA string."""
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100.0


def scan_crispr_guides(seq):
    pam_sites = 0
    usable_guides = []

    # Scan every 3-base window on the forward strand for NGG
    # A 3-base window starting at index i ends at i+3; if i+1 and i+2 are 'GG', it's an NGG PAM.
    for i in range(len(seq) - 2):
        codon_window = seq[i:i + 3]
        if codon_window[1:] == "GG":
            pam_sites += 1
            # We need 20 nucleotides immediately upstream of the PAM (indices i-20 to i)
            if i >= 20:
                guide_seq = seq[i - 20:i]
                gc_content = calculate_gc(guide_seq)
                usable_guides.append({
                    "pam_pos": i + 1,  # 1-based position of the N in NGG
                    "guide": guide_seq,
                    "pam": codon_window,
                    "gc": gc_content,
                    "dist_from_50": abs(gc_content - 50.0)
                })

    # Filter for guides between 40% and 60% GC content (inclusive)
    passed_guides = [g for g in usable_guides if 40.0 <= g["gc"] <= 60.0]

    # Rank by how close they are to 50% GC (then by position as a tie-breaker)
    ranked_guides = sorted(passed_guides, key=lambda x: (x["dist_from_50"], x["pam_pos"]))

    return pam_sites, usable_guides, ranked_guides


# Run the scan
total_pams, usable, passed_ranked = scan_crispr_guides(hbb_sequence)

print("=== SUMMARY METRICS ===")
print(f"Total PAM sites (NGG) on forward strand : {total_pams}")
print(f"Usable guides (full 20bp upstream)      : {len(usable)}")
print(f"Guides passing GC filter (40% - 60%)    : {len(passed_ranked)}\n")

print("=== RANKED GUIDES (40% - 60% GC, closest to 50% first) ===")
print(f"{'Rank':<5} {'PAM Pos':<8} {'20-nt Guide Sequence':<22} {'PAM':<5} {'GC %':<6}")
print("-" * 50)
for idx, g in enumerate(passed_ranked, 1):
    print(f"{idx:<5} {g['pam_pos']:<8} {g['guide']:<22} {g['pam']:<5} {g['gc']:5.1f}%")
