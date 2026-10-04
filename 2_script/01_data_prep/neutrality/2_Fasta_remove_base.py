import argparse
from Bio import SeqIO
import sys

def read_fasta(fasta):
  
    fa_dict = {}
    try:
        with open(fasta, 'r') as fa:
            for seq in SeqIO.parse(fa, "fasta"):
                sequence = str(seq.seq).strip().upper() 
                seq_id = str(seq.id)
                fa_dict[seq_id] = (seq.description, sequence)
        return fa_dict
    except FileNotFoundError:
        print(f"Error: file {fasta} not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)

def remove_codons(sequence, codons_to_remove=None):
   
    if codons_to_remove is None:

        codons_to_remove = ['ATG', 'TGG', 'ATT', 'ATC', 'ATA']

    if len(sequence) % 3 != 0:
        print(f"Warning: sequence length {len(sequence)} is not a multiple of 3; the tail may be truncated")

    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

    filtered_codons = [codon for codon in codons
                      if len(codon) == 3 and codon not in codons_to_remove]

    return ''.join(filtered_codons)

def main():
    parser = argparse.ArgumentParser(
        description='Remove specific triplet codons from CDS sequences. '
                    'By default removes ATG (M), TGG (W) and the isoleucine codons ATT, ATC, ATA.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Codons removed by default:
  ATG - Methionine (Met/M, start codon)
  TGG - Tryptophan (Trp/W)
  ATT - Isoleucine (Ile/I)
  ATC - Isoleucine (Ile/I)
  ATA - Isoleucine (Ile/I)
        """
    )
    parser.add_argument('-a', '--fasta', help='input FASTA file', type=str, required=True)
    parser.add_argument('-o', '--out', help='output file', type=str, required=False)
    parser.add_argument('--codons', help='comma-separated list of codons to remove, e.g.: ATG,TGG,ATT,ATC,ATA',
                       type=str, default='ATG,TGG,ATT,ATC,ATA')

    args = parser.parse_args()

    codons_to_remove = [codon.strip().upper() for codon in args.codons.split(',')]
    print(f"Codons to remove: {codons_to_remove}")

    codon_info = {
        'ATG': 'Methionine (Met/M, start codon)',
        'TGG': 'Tryptophan (Trp/W)',
        'ATT': 'Isoleucine (Ile/I)',
        'ATC': 'Isoleucine (Ile/I)',
        'ATA': 'Isoleucine (Ile/I)'
    }

    print("\nAmino-acid information for removed codons:")
    for codon in codons_to_remove:
        if codon in codon_info:
            print(f"  {codon}: {codon_info[codon]}")

    fasta_dict = read_fasta(args.fasta)
    print(f"\nSuccessfully read {len(fasta_dict)} sequences")

    
    result = {}
    processed_count = 0
    removal_stats = {}  

    for codon in codons_to_remove:
        removal_stats[codon] = 0

    for seq_id, (description, sequence) in fasta_dict.items():
        
        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
        for codon in codons_to_remove:
            removal_stats[codon] += codons.count(codon)

        processed_sequence = remove_codons(sequence, codons_to_remove)
        result[description] = processed_sequence
        processed_count += 1

        
        if processed_count % 100 == 0:
            print(f"Processed {processed_count} sequences...")

    print(f"\nDone. {processed_count} sequences processed")

    print("\nCodon removal statistics:")
    total_removed = 0
    for codon, count in removal_stats.items():
        print(f"  {codon}: {count}")
        total_removed += count
    print(f"Total codons removed: {total_removed}")

    
    if args.out:
        try:
            with open(args.out, "w") as output_file:
                for description, sequence in result.items():
                    output_file.write(f">{description}\n{sequence}\n")
            print(f"\nResults saved to: {args.out}")

            
            original_total_len = sum(len(fasta_dict[seq_id][1]) for seq_id in fasta_dict)
            new_total_len = sum(len(seq) for seq in result.values())
            removed_bases = original_total_len - new_total_len

            print(f"Original total sequence length: {original_total_len} bp")
            print(f"Processed total sequence length: {new_total_len} bp")
            print(f"Bases removed: {removed_bases} bp")
            print(f"Removal proportion: {removed_bases/original_total_len*100:.2f}%")

        except Exception as e:
            print(f"Error writing output file: {e}")
            sys.exit(1)
    else:
       
        for description, sequence in result.items():
            print(f">{description}")
            
            for i in range(0, len(sequence), 80):
                print(sequence[i:i+80])

if __name__ == '__main__':
    main()
