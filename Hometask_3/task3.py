__author__ = 'esenin'

init_array_x = [float(x) / 10 for x in range(10)]
init_array_y = [1.623250, 1.664792, 1.701977, 1.734832, 1.763404, 1.787764, 1.808002, 1.824230, 1.836580, 1.845201]
x_ms = [0.041745, 0.738224, 0.533362]

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


def make_row_out(xys, row, n):
    result = "\t"
    for j in range(n):
        if j % 2 == row % 2 and j <= row <= j + 2 * (len(xys[j]) - 1):
            index = (row - j) / 2
            result += str('{0:0.06f}'.format(xys[j][index]))
            result += "\t\t"
        else:
            result += "\t\t\t"

    return result


def print_out(xys):
    m = len(xys[0])
    n = len(xys)
    for i in range(2 * m):
        print make_row_out(xys, i, n)
    print "\n"


def find_upper_specials(x_m):
    x0 = init_array_y[-1]
    for i in range(1, len(init_array_x)):
        if init_array_y[i] > x_m:
            x0 = init_array_x[i - 1]
            break
    t = (x_m - x0) / (init_array_x[1] - init_array_x[0])
    return x0, t

def find_bottom_specials(x_m):
    x0 = init_array_y[0]
    for i in reversed(range(1, len(init_array_x))):
        if init_array_x[i] < x_m:
            x0 = init_array_x[i - 1]
            break
    t = (x_m - x0) / (init_array_x[1] - init_array_x[0])
    return x0, t


def find_interval(x_m):
    (x0, t_u) = find_upper_specials(x_m)
    (xn, t_d) = find_bottom_specials(x_m)
    print x0, "\t", t_u
    print xn, "\t", t_d


def main():
    result = process_finite_difference(init_array_y)
    print_out(result)
    find_interval(x_ms[0])
    return 0


main()