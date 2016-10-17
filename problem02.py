def sq_hypotenuse(a, b):
    return a ** 2 + b ** 2


if __name__ == '__main__':
    str = input()
    a, b = [int(x) for x in str.split()]
    print(sq_hypotenuse(int(a), int(b)))

