import sys

from Bio import SeqIO

from edit import edit_distance

if __name__ == '__main__':
    records = SeqIO.parse(sys.stdin, format='fasta')
    result = edit_distance(str(r.seq) for r in records)
    print(result.distance, *result.strings, sep='\n')
