from collections import defaultdict

import Bio.Data.CodonTable


if __name__ == "__main__":

    table = Bio.Data.CodonTable.standard_dna_table.forward_table
    back_table = defaultdict(list)
    for k, v in table.items():
        back_table[v].append(k)
    back_table["*"] = Bio.Data.CodonTable.standard_dna_table.stop_codons
    back_table = {k: len(v) for k, v in back_table.items()}

    prot = input()
    prot = prot + "*"  # add missing stop codon

    result = 1
    mod = 1_000_000
    for aa in prot:
        result = result * back_table[aa] % mod

    print(result)
