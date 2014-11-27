# coding=utf-8
__author__ = 'esenin'

import numpy as np


f = lambda x: 1 / (np.exp(x + 0.01) - 0.68)
f_der2 = lambda x: (np.exp(x) * (0.666535 + 0.99005 * np.exp(x))) / (np.exp(x) - 0.673234)**3
bounds = (0.0, 0.4)
a = bounds[0]
b = bounds[1]


def integrate_rectangles(n):
    result = 0.0
    h = float(b - a) / n
    for i in range(1, n + 1):
        result += f(a + float(2*i - 1)/2 * h)
    return result * h


def integrate_trapezoid(n):
    h = float(b - a) / n
    result = (f(a) + f(b)) * 0.5
    for k in range(1, n):
        result += f(a + k * h)
    return result * h


def integrate_by_simpson(n):
    h = float(b - a) / n
    zero_part = (f(a) + f(b))
    first_part = 0.0
    for k in range(1, n):
        first_part += f(a + k * h)
    second_part = 0.0
    for k in range(1, n + 1):
        second_part += f(a + float(2*k - 1)/2 * h)
    return h/6 * (zero_part + 2 * first_part + 4 * second_part)


def exact_value(k, res_n, res_2n):
    return (2**k * res_2n - res_n) / float(2**k - 1)


def accuracy():
    n = 8
    xs = np.linspace(a, b, 100)
    fd2_max = np.max(np.abs(f_der2(xs)))
    return np.abs((b - a)**3) / (24 * n**2) * fd2_max


def process_with_fun(fun, k, name):
    j_8 = fun(8)
    j_16 = fun(16)
    j_r = exact_value(k, j_8, j_16)
    acc = accuracy()
    print_out(j_8, j_16, j_r, acc, name)


def print_out(j1, j2, jr, acc, name):
    format_s = "{:+10.8}"
    print name
    print "\tJ_8 = \t", "\tJ_16 = \t", "\tJ_R = \t", "\tA = \t", "\t|J_8 - J_R| <= A :"
    print format_s.format(j1), \
        format_s.format(j2), \
        format_s.format(jr), \
        format_s.format(acc), \
        "\t\t", np.abs(j1 - jr) <= acc


def main():
    process_with_fun(integrate_rectangles, 2, "Интегрирование по формуле прямоугольников:")
    process_with_fun(integrate_trapezoid, 2, "Интегрирование по формуле трапеции:")
    process_with_fun(integrate_by_simpson, 4, "Интегрирование по формуле Симпсона:")

main()