import sys
from functools import reduce

from Bio import Seq
from Bio import SeqIO

if __name__ == '__main__':
    dna, *introns = SeqIO.parse(sys.stdin, format='fasta')
    s = reduce(lambda acc, intron: acc.replace(str(intron.seq), ''), introns, str(dna.seq))
    print(Seq.translate(s, to_stop=True))
