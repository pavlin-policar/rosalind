import sys
from Bio import SeqIO

N_MATCHING = 3

if __name__ == '__main__':
    records = list(SeqIO.parse(sys.stdin, format='fasta'))
    for record in records:
        for cmp_record in records:
            if record.seq[-N_MATCHING:] == cmp_record.seq[:N_MATCHING] and \
                    record is not cmp_record:
                print(record.id, cmp_record.id)
