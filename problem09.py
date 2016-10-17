if __name__ == '__main__':
    string = input()
    print(string[::-1].translate(str.maketrans('ATCG', 'TAGC')))
