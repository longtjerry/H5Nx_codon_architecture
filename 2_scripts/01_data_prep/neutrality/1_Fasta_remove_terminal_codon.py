import argparse
from Bio import SeqIO


def read_fasta(fasta):
    fa_dict = {}
    with open(fasta,'r') as fa:
        for seq in SeqIO.parse(fa, "fasta"):
            sequence = str(seq.seq).strip()
            seq_id = str(seq.id)
            fa_dict[seq_id] = sequence
    return fa_dict

def cut_faste(fa_dict):
    fasta = {}
    for i in fa_dict.keys():
        sequence = fa_dict.get(i)
        rem = len(sequence) % 3
        if rem == 0:
            seq = sequence[: len(sequence) - 3]
        elif rem == 2:
            seq = sequence[: len(sequence) - 2]
        else:
            seq = sequence[: len(sequence) - 1]
        fasta[i] = seq
    return fasta

def main():
    parser = argparse.ArgumentParser(description='remove terminal codons')
    parser.add_argument('-a', '--fasta', help='input', type=str, required=True)
    parser.add_argument('-o', '--out', help='output', type=str, required=False)
    args = parser.parse_args()

    fasta = read_fasta(args.fasta)
    res = cut_faste(fasta)
    if args.out is None:
        for i in res.keys():
            print(">", i, "\n", res.get(i), sep='')
    else:
        p = open(args.out, "w")
        for i in res.keys():
            p.write(">" + i + "\n" + res.get(i) + "\n")
        p.close()

if __name__ == '__main__':
    main()






