from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from orf import get_orfs

if __name__ == '__main__':
    dna = Seq(input(), alphabet=IUPAC.ambiguous_dna)
    print(max(map(str, get_orfs(dna)), key=len))
