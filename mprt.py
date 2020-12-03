import re
import sys

import requests

if __name__ == "__main__":
    codes = [l.strip() for l in sys.stdin]

    strings = []
    for code in codes:
        res = requests.get("http://www.uniprot.org/uniprot/%s.fasta" % code)
        res = res.content.decode("utf-8")
        split = res.split('\n')
        _, *fasta = split
        fasta = "".join(fasta)
        matches = list(re.finditer(r"(?=N[^P][ST][^P])", fasta))
        if any(matches):
            print(code)
            print(" ".join([str(x.start() + 1) for x in matches]))
