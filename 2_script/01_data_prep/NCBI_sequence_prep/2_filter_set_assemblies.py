import argparse
import os
import sys

SEGMENTS = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']


def filter_file(inpath, outpath):
    kept = 0
    removed = 0
    with open(inpath, 'r') as fin, open(outpath, 'w') as fout:
        header = None
        seq_lines = []
        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    assembly = header.split('|')[-1].strip()
                    if assembly.startswith('set:'):
                        removed += 1
                    else:
                        fout.write(header + '\n')
                        fout.write('\n'.join(seq_lines) + '\n')
                        kept += 1
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            assembly = header.split('|')[-1].strip()
            if assembly.startswith('set:'):
                removed += 1
            else:
                fout.write(header + '\n')
                fout.write('\n'.join(seq_lines) + '\n')
                kept += 1
    return kept, removed


def main():
    parser = argparse.ArgumentParser(
        description="Remove 'set:' assembly entries from segment FASTA files"
    )
    parser.add_argument('--indir', '-i', required=True, help='Input segment directory')
    parser.add_argument('--outdir', '-o', help='Output directory (omit for --inplace)')
    parser.add_argument('--inplace', action='store_true', help='Overwrite files in-place')
    args = parser.parse_args()

    if not args.inplace and not args.outdir:
        print("Error: specify --outdir or --inplace")
        sys.exit(1)

    outdir = args.indir if args.inplace else args.outdir
    os.makedirs(outdir, exist_ok=True)

    total_kept = 0
    total_removed = 0

    print(f"{'Segment':<10} {'Kept':>8} {'Removed':>8}")
    print("-" * 28)

    for seg in SEGMENTS:
        infile = os.path.join(args.indir, f"{seg}.fasta")
        if not os.path.exists(infile):
            continue
        outfile = os.path.join(outdir, f"{seg}.fasta")
        kept, removed = filter_file(infile, outfile)
        total_kept += kept
        total_removed += removed
        print(f"{seg:<10} {kept:>8} {removed:>8}")

    print("-" * 28)
    print(f"{'TOTAL':<10} {total_kept:>8} {total_removed:>8}")
    print(f"\nCleaned files written to: {outdir}/")


if __name__ == '__main__':
    main()
