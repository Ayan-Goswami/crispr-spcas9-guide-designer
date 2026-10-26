### CRISPR SpCas9 Guide RNA Designer - Standalone Reproducible Script

SEQUENCE = (
    "ATGGAGGGCCATGTCAAGCGCCCCATGAATGCATTTATGGTGTGGTCCCGTGGTGAGAGGCACAAGTTGG"
    "CCCAGCAGAATCCCAGCATGCAAAATACAGAGATCAGCAAGCAGCTGGGATGCAGGTGGAAAAGCCTTAC"
    "AGAAGCCGAAAAAAGGCCCTTTTTCCAGGAGGCACAGAGATTGAAGATCCTACACAGAGAGAAATACCCA"
    "AACTATAAATATCAGCCTCATCGGAGGGCTAAAGTGTCACAGAGGAGTGGCATTTTACAGCCTGCAGTTG"
    "CCTCAACAAAACTGTACAACCTTCTGCAGTGGGACAGGAACCCACATGCCATCACATACAGGCAAGACTG"
    "GAGTAGAGCTGCACACCTGTACTCCAAAAACCAGCAAAGCTTTTATTGGCAGCCTGTTGATATCCCCACT"
    "GGGCACCTGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTCCATAACCACCACCAGCAGCAACAGC"
    "AGTTCTATGACCACCACCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTCCATGACCACCACCA"
    "GCAGAAGCAGCAGTTTCATGACCACCACCAGCAGCAACAGCAGTTCCATGACCACCACCACCACCACCAG"
    "GAGCAGCAGTTCCATGACCACCACCAGCAGCAACAGCAGTTCCATGACCACCAGCAGCAGCAGCAGCAGC"
    "AGCAGCAGCAGCAGTTCCATGACCACCACCAGCAGAAGCAGCAGTTCCATGACCACCACCACCACCAACA"
    "GCAGCAGCAGTTCCATGACCACCAGCAGCAGCAGCAGCAGTTCCATGACCACCAGCAGCAGCAGCATCAG"
    "TTCCATGACCACCCCCAGCAGAAGCAGCAGTTCCATGACCACCCCCAGCAGCAACAGCAGTTCCATGACC"
    "ACCACCACCAGCAGCAGCAGAAGCAGCAGTTCCATGACCACCACCAGCAGAAGCAGCAGTTCCATGACCA"
    "CCACCAGCAGAAGCAGCAGTTCCATGACCACCACCAGCAGCAACAGCAGTTCCATGACCACCACCAGCAG"
    "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTCCACGACCAGCAGCTTACCTACTTACTAACAGCTG"
    "ACATCACTGGTGAGCATACACCATACCAGGAGCACCTCAGCACAGCCCTGTGGTTGGCAGTCTCATGA"
)

STRAND_MODE = "+"   # Options: "both", "+", "-"
MIN_GC = 40
MAX_GC = 60
TARGET_POS = 20       # 1-based target coordinate (or None)
SORT_BY = "gc"           # Options: "gc", "target", "position"

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp[b] for b in reversed(seq))

def calculate_gc(seq):
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100.0

def scan_crispr_guides(seq, strand_mode='both', min_gc=40.0, max_gc=60.0, target_pos=None, sort_by='gc'):
    seq_len = len(seq)
    usable_guides = []
    fwd_pams = rev_pams = 0

    # 1. Forward strand (+) scan (NGG)
    if strand_mode in ('both', '+'):
        for i in range(seq_len - 2):
            if seq[i+1:i+3] == 'GG':
                fwd_pams += 1
                if i >= 20:
                    guide = seq[i-20:i]
                    gc = calculate_gc(guide)
                    cut_site = i - 2  # 1-based: 3 bp upstream of PAM
                    usable_guides.append({
                        'strand': '+',
                        'guide_start': i - 19,
                        'guide_end': i,
                        'pam_pos': i + 1,
                        'cut_site': cut_site,
                        'dist_target': abs(cut_site - target_pos) if target_pos is not None else None,
                        'guide': guide,
                        'pam': seq[i:i+3],
                        'gc': gc,
                        'dist_50': abs(gc - 50.0)
                    })

    # 2. Reverse strand (-) scan (CCN on forward strand)
    if strand_mode in ('both', '-'):
        for i in range(seq_len - 2):
            if seq[i:i+2] == 'CC':
                rev_pams += 1
                if i + 23 <= seq_len:
                    guide = reverse_complement(seq[i+3:i+23])
                    gc = calculate_gc(guide)
                    cut_site = i + 6  # 1-based: 3 bp inside protospacer from PAM
                    usable_guides.append({
                        'strand': '-',
                        'guide_start': i + 4,
                        'guide_end': i + 23,
                        'pam_pos': i + 1,
                        'cut_site': cut_site,
                        'dist_target': abs(cut_site - target_pos) if target_pos is not None else None,
                        'guide': guide,
                        'pam': reverse_complement(seq[i:i+3]),
                        'gc': gc,
                        'dist_50': abs(gc - 50.0)
                    })

    passed = [g for g in usable_guides if min_gc <= g['gc'] <= max_gc]

    if sort_by == 'target' and target_pos is not None:
        passed.sort(key=lambda x: (x['dist_target'], x['dist_50'], x['pam_pos']))
    elif sort_by == 'position':
        passed.sort(key=lambda x: (x['pam_pos'], x['dist_50']))
    else:
        passed.sort(key=lambda x: (x['dist_50'], x['pam_pos'], 0 if x['strand'] == '+' else 1))

    return fwd_pams + rev_pams, usable_guides, passed

total_pams, usable, passed = scan_crispr_guides(
    SEQUENCE, STRAND_MODE, MIN_GC, MAX_GC, TARGET_POS, SORT_BY
)

print(f'Total PAMs ({STRAND_MODE}): {total_pams} | Usable: {len(usable)} | Passed ({MIN_GC}-{MAX_GC}% GC): {len(passed)}')
print('-' * 82)
cols = ['Rank', 'Str', 'Span', 'PAM Pos', 'Cut', 'dTarget', '20-nt Guide', 'PAM', 'GC%']
print(f'{cols[0]:<5} {cols[1]:<4} {cols[2]:<10} {cols[3]:<8} {cols[4]:<6} {cols[5]:<8} {cols[6]:<22} {cols[7]:<5} {cols[8]:<6}')
print('-' * 82)
for idx, g in enumerate(passed, 1):
    span = f"{g['guide_start']}-{g['guide_end']}"
    dt = f"{g['dist_target']}bp" if g['dist_target'] is not None else 'N/A'
    row = (
        f"{idx:<5} "
        f"{g['strand']:<4} "
        f"{span:<10} "
        f"{g['pam_pos']:<8} "
        f"{g['cut_site']:<6} "
        f"{dt:<8} "
        f"{g['guide']:<22} "
        f"{g['pam']:<5} "
        f"{g['gc']:5.1f}%"
    )
    print(row)
