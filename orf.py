import sys
from re import findall

from Bio import SeqIO
from Bio.Seq import translate


def get_orfs(seq):
    orfs = set()
    for strand in (seq, seq.reverse_complement()):
        for offset in range(3):
            matches = findall(
                r'(?=(ATG(?:(?:\w{3})*?)(?:TAG|TGA|TAA)))',
                str(strand)[offset:])
            for match in matches:
                orfs.add(translate(match, to_stop=True))
    return orfs


if __name__ == '__main__':
    record, = list(SeqIO.parse(sys.stdin, 'fasta'))
    print(*get_orfs(record.seq), sep='\n')
