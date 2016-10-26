from itertools import chain

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from re import finditer

if __name__ == '__main__':
    dna = Seq(input(), alphabet=IUPAC.ambiguous_dna)

    orfs = [[strand[x.start():].translate(to_stop=True)
             for x in finditer('ATG', str(strand))]
            for strand in (dna, dna.reverse_complement())]
    orfs = {str(orf) for orf in chain(*orfs)}

    print(max(map(str, orfs), key=len))
