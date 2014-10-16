__author__ = 'esenin'


init_array = [1.623250, 1.664792, 1.701977, 1.734832, 1.763404, 1.787764, 1.808002, 1.824230, 1.836580, 1.845201]
EPS = 10**(-5)


def calc_finite_difference(xs):
    result = []
    for i in range(len(xs) - 1):
        result.append(xs[i] - xs[i + 1])
    return result


def has_nonzero_values(xs):
    return any([abs(i) > EPS for i in xs])


def process_finite_difference(xs):
    result = [xs]
    current = xs
    counter = 0
    while has_nonzero_values(current):
        current = calc_finite_difference(current)
        result.append(current)
        counter += 1
        assert counter < 1000, "\tIterations count exceeded.Stop"

    return result


def make_row_out(xys, row, m, n):
    result = "\t"
    for j in range(n):
        if j % 2 == row % 2 and j <= row <= j + 2 * (len(xys[j]) - 1):
            index = (row - j) / 2
            #result += str(xys[j][(row - j) / 2])
            result += str('{0:0.05f}'.format(xys[j][index])
            result += "\t\t"
        else:
            result += "\t\t\t"

    return result


def print_out(xys):
    m = len(xys[0])
    n = len(xys)
    for i in range(2 * m):
        print make_row_out(xys, i, m, n)


def main():
    result = process_finite_difference(init_array)
    print_out(result)
    return 0


main()