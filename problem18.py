import sys
from itertools import chain
from re import finditer

from Bio import SeqIO

if __name__ == '__main__':
    record, = list(SeqIO.parse(sys.stdin, 'fasta'))
    dna = record.seq

    orfs = [[strand[x.start():].translate(to_stop=True)
             for x in finditer('ATG', str(strand))]
            for strand in (dna, dna.reverse_complement())]
    orfs = {str(orf) for orf in chain(*orfs)}

    print(*orfs, sep='\n')
