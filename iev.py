if __name__ == "__main__":
    s = input()
    vals = list(map(int, s.split(" ")))
    exp = vals[0] + vals[1] + vals[2] + 0.75 * vals[3] + 0.5 * vals[4]
    print(exp * 2)
