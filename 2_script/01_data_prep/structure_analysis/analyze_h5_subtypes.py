import os
import pandas as pd
from collections import defaultdict

def read_fasta(filename):
   
    sequences = []
    current_seq = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_seq:
            sequences.append(''.join(current_seq))
    return sequences

def analyze_subtype_files(file_list):
   
    subtype_sequences = {}
    all_lengths = []

    for fname, subtype in file_list:
        if not os.path.exists(fname):
            print(f"Warning: {fname} not found. Skipping.")
            continue
        seqs = read_fasta(fname)
        if not seqs:
            print(f"Warning: No sequences in {fname}")
            continue
        subtype_sequences[subtype] = seqs
        
        for s in seqs:
            all_lengths.append(len(s))
        print(f"{subtype}: {len(seqs)} sequences")

    if not subtype_sequences:
        print("No data loaded. Exiting.")
        return

    unique_lengths = set(all_lengths)
    if len(unique_lengths) > 1:
        print(f"Warning: Sequences have different lengths: {unique_lengths}")
        print("Using the most common length (assuming alignment gaps).")
       
        from collections import Counter
        aln_length = Counter(all_lengths).most_common(1)[0][0]
    else:
        aln_length = list(unique_lengths)[0]

    if aln_length % 3 != 0:
        print(f"Warning: Alignment length {aln_length} not multiple of 3.")
        print("Make sure sequences are codon-aligned (CDS only).")

    print(f"\nAlignment length: {aln_length} nt -> {aln_length//3} amino acids")

    arg_codons = ['AGA', 'AGG', 'CGT', 'CGC', 'CGA', 'CGG']

    results = []

    for codon_start in range(0, aln_length, 3):
        aa_pos = codon_start // 3 + 1

        for subtype, seqs in subtype_sequences.items():
          
            counts = {codon: 0 for codon in arg_codons}
            total_valid = 0   
            total_arg = 0

            for seq in seqs:
                if codon_start + 3 > len(seq):
                    continue
                codon = seq[codon_start:codon_start+3]
                if '-' in codon or 'N' in codon:
                    continue
                total_valid += 1
                if codon in arg_codons:
                    counts[codon] += 1
                    total_arg += 1

            if total_arg > 0:
                aga_agg_ratio = (counts['AGA'] + counts['AGG']) / total_arg
                aga_prop = counts['AGA'] / total_arg
                agg_prop = counts['AGG'] / total_arg
            else:
                aga_agg_ratio = 0.0
                aga_prop = 0.0
                agg_prop = 0.0

            coverage = total_valid / len(seqs) if len(seqs) > 0 else 0


            if total_arg > 0:
                results.append({
                    'AA_Position': aa_pos,
                    'Subtype': subtype,
                    'AGA_Count': counts['AGA'],
                    'AGG_Count': counts['AGG'],
                    'CGT_Count': counts['CGT'],
                    'CGC_Count': counts['CGC'],
                    'CGA_Count': counts['CGA'],
                    'CGG_Count': counts['CGG'],
                    'Total_Arg_Codons': total_arg,
                    'Total_Valid_Codons': total_valid,
                    'Subtype_Sequences': len(seqs),
                    'Coverage': coverage,
                    'AGA_AGG_Ratio': aga_agg_ratio,
                    'AGA_Proportion': aga_prop,
                    'AGG_Proportion': agg_prop
                })

    df = pd.DataFrame(results)
    if df.empty:
        print("No arginine positions found.")
        return

    df = df.sort_values(['AA_Position', 'Subtype'])
    out_csv = 'h5_ha_arg_by_subtype.csv'
    df.to_csv(out_csv, index=False)
    print(f"\nResults saved to {out_csv}")
    print(f"Total rows: {len(df)}")
    return df

if __name__ == "__main__":
    files = [
        ('ORF_HA_H5N1.fasta', 'H5N1'),
        ('ORF_HA_H5N2.fasta', 'H5N2'),
        ('ORF_HA_H5N6.fasta', 'H5N6'),
        ('ORF_HA_H5N8.fasta', 'H5N8'),
        ('ORF_HA_Others.fasta', 'Others')
    ]
    analyze_subtype_files(files)
