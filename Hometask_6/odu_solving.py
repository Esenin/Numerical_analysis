# coding=utf-8
__author__ = 'esenin'

'''
    Koshi problem solution
    {y' = 1 - sin(ax + y) + ky / (2+x)
    {y(0) = 0   [0, 1]
        a = 1 (0.25) 2,  k = -0.3 (0.2) 0.5,  h = 0.1

'''

import numpy as np
import matplotlib.pyplot as plt

# ######################################################################################################################
#
#   CONSTANTS
#
#######################################################################################################################
format_s = "{:+12.8}"
startPoint = (0.0, 0.0)
bounds = (0.0, 1.0)
fun = lambda a, k, x, yx: 1 - np.sin(a * x + yx) + k * yx / (2 + x)
a_koef = np.linspace(1.0, 2.0, 4)
k_koef = np.linspace(-0.3, 0.5, 5)
h_general = 0.1

EPS = 10 ** (-5)

fun_der_x = lambda a, k, x, y: -a * np.cos(a*x + y) - k * y / ((x + 2)**2)
fun_der_y = lambda a, k, x, y: k / (x + 2) - np.cos(a*x + y)


#######################################################################################################################
#
#   HELPERS
#
#######################################################################################################################
def calc_xarray(h):
    xs = np.linspace(bounds[0], bounds[1], np.round(bounds[1] - bounds[0]) / h)
    assert xs[0] == startPoint[0], "xs must be started with left bound value"
    return xs


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


def print_out_finite_difference(xys):
    print " x\t\t y"
    m = len(xys[0])
    n = len(xys)
    for i in range(2 * m):
        print make_row_out(xys, i, n)
    print "\n"


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
        assert counter < 100, "\tIterations count exceeded. Stop"

    return result


def calc_finite_difference(xs):
    result = []
    for i in range(len(xs) - 1):
        result.append(xs[i + 1] - xs[i])
    return result


#######################################################################################################################
#
#   ODU SOLUTIONS
#
#######################################################################################################################
def euler_solution(a, k, h):
    xs = calc_xarray(h)
    ys = [startPoint[1]]
    for x in xs:
        ys.append(ys[-1] + h * fun(a, k, x, ys[-1]))
    ys = ys[0: -1]
    assert len(xs) == len(ys), "xs, ys size problem"
    return xs, ys


def rungekutta_solution(a, k, h=h_general):
    xs = calc_xarray(h)
    ys = [startPoint[1]]
    for x in xs:
        yk = ys[-1]
        k1 = h * fun(a, k, x, yk)
        k2 = h * fun(a, k, x + h / 2.0, yk + k1 / 2)
        k3 = h * fun(a, k, x + h / 2.0, yk + k2 / 2)
        k4 = h * fun(a, k, x + h, yk + k3)
        ys.append(yk + 1.0 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4))
    ys = ys[0: -1]
    assert len(xs) == len(ys), "xs, ys size problem"
    return xs, ys


def adams_solution(a, k, h=h_general):
    xys_r = rungekutta_solution(a, k, h)
    xs = calc_xarray(h)
    ys = xys_r[1][0:4]
    ettas = []
    for i, x in enumerate(xys_r[0]):
        ettas.append(h * fun(a, k, x, xys_r[1][i]))
    fin_dif = process_finite_difference(ettas)
    for i in range(3, len(xs)):
        yk = ys[-1]
        ys.append(yk + ettas[i] + 0.5 * fin_dif[1][i - 1] + 5.0 / 12.0 * fin_dif[2][i - 2] + 3.0 / 8.0 * fin_dif[3][
            i - 3] + 251.0 / 720.0 * fin_dif[4][i - 4])
    ys = ys[0: -1]
    assert len(xs) == len(ys), "xs, ys size problem"
    return xs, ys


#######################################################################################################################
#
#   GENERAL
#
#######################################################################################################################
def calc_accuracy(a, k, adams_res, iterator):
    acc = 400
    xs = np.linspace(bounds[0], bounds[1], acc)
    ys = np.linspace(np.min(adams_res[1]), np.max(adams_res[1]), acc)
    m1 = np.max(np.abs(fun(a, k, xs, ys)))
    m2 = np.max(np.abs(fun_der_x(a, k, xs, ys)))
    m3 = np.max(np.abs(fun_der_y(a, k, xs, ys)))
    m4 = m2 + m1 * m3
    return (m4 * h_general / 2 * m3) * np.exp(m3 * (adams_res[0][0] - adams_res[0][iterator]))


def process_euler(a, k, h):
    adams_res = adams_solution(a, k)
    eulers_res = euler_solution(a, k, h)
    print "\n----------------------------------------------------------------------------"
    print "\t\tМетод Эйлера\n\ta = ", a, "\tk = ", k, "\n"
    print "\tx_k\t\t\t\ty_k\t\t\t|Y_Adams - Y_Euler| \t<=\t\t\tacc"
    for i, x in enumerate(eulers_res[0]):
        xk = eulers_res[0][i]
        y_e = eulers_res[1][i]
        y_a = adams_res[1][i]
        acc = calc_accuracy(a, k, adams_res, i)
        delta = np.abs(y_a - y_e)
        print format_s.format(xk), \
            format_s.format(y_e), "\t\t\t", \
            format_s.format(delta), "\t\t", \
            delta <= acc, "\t", \
            format_s.format(acc)
    print "\n----------------------------------------------------------------------------"
    return eulers_res


def process_runge_and_adams(a, k):
    runge_res = rungekutta_solution(a, k)
    adams_res = adams_solution(a, k)
    print "\n----------------------------------------------------------------------------"
    print "\t\tМетод Рунге-Кутта  \t\t Метод Адамса\n\ta = ", a, "\tk = ", k, "\n"
    print "\tx_k\t\t\t\ty_k(Runge)\t\t\ty_k(Adams)"
    for i, x in enumerate(runge_res[0]):
        xk = runge_res[0][i]
        print format_s.format(xk), "\t\t", \
            format_s.format(runge_res[1][i]), "\t\t", \
            format_s.format(adams_res[1][i])
    print "\nk=10   |Y(Runge) - Y(Adams)|\t<=\t\t\tacc"
    acc = calc_accuracy(a, k, adams_res, len(adams_res[0]) - 1)
    delta = np.abs(adams_res[1][-1] - runge_res[1][-1])
    print "\t\t", format_s.format(delta), "\t\t\t", delta <= acc, "\t\t", acc
    print "\n----------------------------------------------------------------------------"
    return runge_res, adams_res


def main():
    euler_res = process_euler(a_koef[0], k_koef[0], h_general)
    eunge_res, adams_res = process_runge_and_adams(a_koef[0], k_koef[0])


    adams_line, = plt.plot(adams_res[0], adams_res[1], label='Adams')
    eulers_line, = plt.plot(euler_res[0], euler_res[1], label='Euler')
    plt.legend()
    plt.show()


main()
