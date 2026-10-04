import argparse
import os
import re
import sys
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:
    print("Error: biopython required. Install: pip install biopython")
    sys.exit(1)


SEGMENT_PATTERNS = {
    'PB2': [r'\bPB2\b', r'polymerase basic 2'],
    'PB1': [r'\bPB1\b(?!-F2)', r'polymerase basic 1'],
    'PA':  [r'\bPA\b(?!-X)', r'polymerase acidic'],
    'HA':  [r'\bHA\b', r'hemagglutinin', r'Hemagglutinin'],
    'NP':  [r'\bNP\b', r'nucleoprotein', r'Nucleoprotein'],
    'NA':  [r'\bNA\b', r'neuraminidase', r'Neuraminidase'],
    'MP':   [r'\bmatrix\b', r'Matrix', r'\bM1\b', r'\bM2\b'],
    'NS':  [r'\bNS\b', r'non.structural', r'\bNEP\b', r'\bNS1\b'],
}

SEGMENT_ORDER = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']



def extract_isolate(header):
   
    idx = header.find('(A/')
    if idx == -1:
        m = re.search(r'\(([A-Z]/[^)]+)\)', header)
        return m.group(1) if m else None

    idx += 1  
    depth = 1
    for i in range(idx + 2, len(header)):
        if header[i] == '(':
            depth += 1
        elif header[i] == ')':
            depth -= 1
            if depth == 0:
                return header[idx:i]
    return None


def extract_assembly(header):
   
    m = re.search(r'\|(\S+)', header)
    return m.group(1) if m else None


def identify_segment(header):

    for seg in SEGMENT_ORDER:
        for pat in SEGMENT_PATTERNS[seg]:
            if re.search(pat, header, re.IGNORECASE):
                return seg
    return 'UNKNOWN'


def process(fasta_path, outdir):
    os.makedirs(outdir, exist_ok=True)

    handles = {}
    for seg in SEGMENT_ORDER + ['UNKNOWN']:
        fpath = os.path.join(outdir, f"{seg}.fasta")
        handles[seg] = open(fpath, 'w')

    stats = defaultdict(int)
    discarded = 0

    print(f"Processing: {fasta_path}")
    print("-" * 50)

    for record in SeqIO.parse(fasta_path, 'fasta'):
        header = record.description
        seq = str(record.seq)

        assembly = extract_assembly(header)
        if not assembly:
            discarded += 1
            continue

        isolate = extract_isolate(header)
        if not isolate:
       
            isolate = header.split()[0] if header else 'unknown'

        seg = identify_segment(header)
        stats[seg] += 1

        new_header = f">{isolate}|{assembly}".replace(' ', '_')
        handles[seg].write(f"{new_header}\n{seq}\n")

    for h in handles.values():
        h.close()

    print(f"\n{'Segment':<10} {'Count':>8}")
    print("-" * 20)
    total_kept = 0
    for seg in SEGMENT_ORDER + ['UNKNOWN']:
        cnt = stats[seg]
        total_kept += cnt
        marker = " ✓" if seg in SEGMENT_ORDER and cnt > 0 else ""
        print(f"{seg:<10} {cnt:>8}{marker}")
    print("-" * 20)
    print(f"{'TOTAL':<10} {total_kept:>8}")
    print(f"\nDiscarded (no assembly): {discarded}")
    print(f"Output directory: {outdir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Split meta influenza FASTA into 8 segment files."
    )
    parser.add_argument('--fasta', '-f', required=True, help='Input meta FASTA')
    parser.add_argument('--outdir', '-o', default='./segments', help='Output dir')
    args = parser.parse_args()
    process(args.fasta, args.outdir)


if __name__ == '__main__':
    main()
