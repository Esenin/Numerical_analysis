# coding=utf-8
__author__ = 'esenin'

import math

array_x = [float(x) / 10 for x in range(10)]
array_y = [1.623250, 1.664792, 1.701977, 1.734832, 1.763404, 1.787764, 1.808002, 1.824230, 1.836580, 1.845201]
x_ms = [0.041745, 0.738224, 0.533362]
y_ms = 1.766753

EPS = 10**(-5)


def calc_finite_difference(xs):
    result = []
    for i in range(len(xs) - 1):
        result.append(xs[i + 1] - xs[i])
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
        assert counter < 100, "\tIterations count exceeded.Stop"

    return result


def make_row_out(xys, row, n):
    result = ""
    if row % 2 == 0:
        result = str(array_x[row / 2])
    else:
        result += "\t"
    result += " \t"


    for j in range(n):
        if j % 2 == row % 2 and j <= row <= j + 2 * (len(xys[j]) - 1):
            index = (row - j) / 2
            result += str('{0:0.06f}'.format(xys[j][index]))
            result += "\t"
        else:
            result += "\t\t\t"

    return result


def print_out(xys):
    print " x\t\t y"
    m = len(xys[0])
    n = len(xys)
    for i in range(2 * m):
        print make_row_out(xys, i, n)
    print "\n"


def polynom_at_begin(xys):
    t = (x_ms[0] - array_x[0]) / (array_x[1] - array_x[0])
    factor = 1
    t_op = t
    result = array_y[0]
    for i in range(1, 5):
        result += (1 / math.factorial(factor)) * xys[i][0] * t_op
        factor += 1
        t_op *= t_op - 1

    return t, result


def polynom_at_end(xys):
    n_id = 8
    t = (x_ms[1] - array_x[n_id]) / (array_x[1] - array_x[0])
    factor = 1
    t_op = t
    result = array_y[n_id]
    for i in range(1, 5):
        n_id -= 1
        result += (1 / math.factorial(factor)) * xys[i][n_id] * t_op
        factor += 1
        t_op *= t_op + 1

    return t, result


def polynom_at_middle(xys):
    n_id = 5
    t = (x_ms[2] - array_x[n_id]) / (array_x[1] - array_x[0])
    factor = 1
    t_op = t
    result = array_y[n_id]
    for i in range(1, 5):
        n_id -= 1 if i % 2 == 0 else 0
        result += (1 / math.factorial(factor)) * xys[i][n_id] * t_op
        factor += 1
        t_op *= (t_op - 1) if i % 2 == 1 else (t_op + 1)

    return t, result


def reverse_polynom(xys):
    print "Обратное интерполирование\n"
    print " i\t\t\tt\t\t\t\tu(t)"
    n_id = 4
    t_k = 0
    t_kn = 1
    counter = 0
    while abs(t_k - t_kn) > EPS:
        counter += 1
        t_kn = y_ms
        t_kn -= array_y[n_id]
        t_kn -= t_k * (t_k - 1) / 2 * xys[2][n_id]
        t_kn -= t_k * (t_k - 1) * (t_k - 2) / 6 * xys[3][n_id]
        t_kn -= t_k * (t_k - 1) * (t_k - 2) * (t_k - 3) / math.factorial(4) * xys[4][n_id]
        t_kn *= 1 / (xys[1][n_id])

        print counter, "\t\t", "{:0.10f}".format(t_k), "\t\t", "{:10.10f}".format(t_kn)
        (t_k, t_kn) = (t_kn, t_k)
    print "\nx* = ", t_k * (array_x[1] - array_x[0]) + array_x[n_id]


def main():
    result = process_finite_difference(array_y)
    print_out(result)

    (t1, pt1) = polynom_at_begin(result)
    print "Интерполирование в начале/конце/середине таблицы:"
    print "x*_1 = ", x_ms[0], "\tt = ", t1, "\tP(t) = ", pt1
    (t2, pt2) = polynom_at_end(result)
    print "x*_2 = ", x_ms[1], "\tt = ", t2, "\tP(t) = ", pt2
    (t3, pt3) = polynom_at_middle(result)
    print "x*_3 = ", x_ms[2], "\tt = ", t3, "\tP(t) = ", pt3, "\n"

    reverse_polynom(result)

    return 0


main()