import argparse
import sys
import os
import glob


def rna_to_dna(seq):
    return seq.replace('U', 'T').replace('u', 't')


def process_fasta(in_path, out_path):
    line_count = 0
    seq_count = 0
    u_count = 0

    with open(in_path, 'r', encoding='utf-8') as fin,          open(out_path, 'w', encoding='utf-8') as fout:

        for line in fin:
            line_count += 1
            stripped = line.rstrip('\n')

            if stripped.startswith('>'):
                fout.write(stripped + '\n')
                seq_count += 1
            else:
                converted = rna_to_dna(stripped)
                u_count += converted.count('T') - stripped.count('T')
                u_count += stripped.count('U') + stripped.count('u')
                fout.write(converted + '\n')

    return line_count, seq_count, u_count


def main():
    parser = argparse.ArgumentParser(
        description='Convert RNA sequences to DNA in FASTA file(s).'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA file (or glob pattern for batch)')
    parser.add_argument('-o', '--output',
                        help='Output FASTA file (required for single file, unless --in-place)')
    parser.add_argument('-d', '--outdir',
                        help='Output directory (for batch mode)')
    parser.add_argument('--in-place', action='store_true',
                        help='Overwrite original file(s)')

    args = parser.parse_args()

    if '*' in args.input or '?' in args.input:
        files = sorted(glob.glob(args.input))
        if not files:
            print(f"Error: No files matched pattern: {args.input}")
            sys.exit(1)
    else:
        files = [args.input]

    if len(files) == 1 and not args.in_place and not args.output:
        print("Error: Single file mode requires -o/--output or --in-place")
        sys.exit(1)

    if len(files) > 1 and not args.in_place and not args.outdir:
        print("Error: Batch mode requires -d/--outdir or --in-place")
        sys.exit(1)

    total_files = 0
    for filepath in files:
        if not os.path.exists(filepath):
            print(f"Warning: File not found, skipping: {filepath}")
            continue

        if args.in_place:
            out_path = filepath + '.tmp'
        elif len(files) == 1:
            out_path = args.output
        else:
            os.makedirs(args.outdir, exist_ok=True)
            basename = os.path.basename(filepath)
            out_path = os.path.join(args.outdir, basename)

        lines, seqs, u_replaced = process_fasta(filepath, out_path)

        if args.in_place:
            os.replace(out_path, filepath)
            out_path = filepath

        total_files += 1
        print(f"[{total_files}] {filepath} -> {out_path}")
        print(f"    Sequences: {seqs}, Lines: {lines}, U/u replaced: {u_replaced}")

    print(f"\nDone. {total_files} file(s) processed.")


if __name__ == '__main__':
    main()
