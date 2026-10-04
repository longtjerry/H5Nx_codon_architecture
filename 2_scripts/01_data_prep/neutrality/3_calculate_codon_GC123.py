import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction as GC

def calculate_GC123(seq):
    """Calculate GC content at the first, second, and third codon positions"""
    GC1 = []
    GC2 = []
    GC3 = []

    seq_len = len(seq)
    if seq_len % 3 != 0:
        print(f"Warning: sequence length {seq_len} is not a multiple of 3; truncating to {seq_len - (seq_len % 3)}")
        seq = seq[:seq_len - (seq_len % 3)]

    for i in range(0, len(seq), 3):
        if i + 2 < len(seq):
            GC1.append(str(seq[i]))      
            GC2.append(str(seq[i+1]))    
            GC3.append(str(seq[i+2]))

    GC1_content = round(GC(''.join(GC1)) * 100, 2) if GC1 else 0
    GC2_content = round(GC(''.join(GC2)) * 100, 2) if GC2 else 0
    GC3_content = round(GC(''.join(GC3)) * 100, 2) if GC3 else 0

    return GC1_content, GC2_content, GC3_content

def main():
 
    if len(sys.argv) != 3:
        print("Usage: python calculate_codon_GC123.py cds.fasta output.txt")
        print("Arguments:")
        print("  cds.fasta   : input FASTA file of CDS sequences")
        print("  output.txt  : output results file")
        sys.exit(1)

    in_fasta = sys.argv[1]
    out_results = sys.argv[2]

    try:
        with open(in_fasta, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: input file {in_fasta} not found")
        sys.exit(1)

    results = []
    total_sequences = 0
    sequences_with_issues = 0

    for rec in SeqIO.parse(in_fasta, 'fasta'):
        total_sequences += 1
        seq = rec.seq.upper()  

        valid_bases = set('ATCGN')
        if not all(base in valid_bases for base in seq):
            print(f"Warning: sequence {rec.id} contains invalid characters; skipped")
            sequences_with_issues += 1
            continue

        try:
            GC_total = round(GC(seq) * 100, 2)
            GC1s, GC2s, GC3s = calculate_GC123(seq)
            GC12 = round((GC1s + GC2s) / 2, 2)

            results.append({
                'id': rec.id,
                'GC_total': GC_total,
                'GC1s': GC1s,
                'GC2s': GC2s,
                'GC3s': GC3s,
                'GC12': GC12
            })
        except Exception as e:
            print(f"Error processing sequence {rec.id}: {e}")
            sequences_with_issues += 1
            continue

    try:
        with open(out_results, 'w') as fw:
           
            fw.write("gene_name\tGC_total\tGC1s\tGC2s\tGC3s\tGC12\n")

            for result in results:
                fw.write(f"{result['id']}\t{result['GC_total']}\t{result['GC1s']}\t{result['GC2s']}\t{result['GC3s']}\t{result['GC12']}\n")

        print(f"Done!")
        print(f"Total sequences: {total_sequences}")
        print(f"Successfully processed: {len(results)}")
        print(f"Problematic sequences: {sequences_with_issues}")
        print(f"Results saved to: {out_results}")

    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
