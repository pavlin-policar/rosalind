import re

import requests

if __name__ == "__main__":
    accession = input()

    res = requests.get(f"http://www.uniprot.org/uniprot/{accession}.txt")
    res = res.content.decode("utf-8")

    matches = re.findall(r"GO; GO:.*:(.+);", res)
    print(*matches, sep="\n")
