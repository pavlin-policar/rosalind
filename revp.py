import sys

from Bio import SeqIO

if __name__ == '__main__':
    record, = list(SeqIO.parse(sys.stdin, format='fasta'))
    s = record.seq
    for l in range(4, 13):
        for i in range(len(s) - l + 1):
            if s[i:i + l] == s[i:i + l].reverse_complement():
                print(i + 1, l)
