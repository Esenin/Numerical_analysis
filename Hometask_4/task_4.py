# coding=utf-8
__author__ = 'esenin'

''' Lagrange polynomial '''

import numpy as np


f = lambda xs: (np.cos(2 * xs))**2
g = lambda xs: np.sin(8 * xs)

f_derivate_7 = lambda x: 8192 * np.sin(4 * x)
f_derivate_13 = lambda x: -33554432 * np.sin(4 * x)
g_derivate_7 = lambda x: -2097152 * np.cos(8 * x)
g_derivate_13 = lambda x: 549755813888 * np.cos(8 * x)


left_bound = 0
right_bound = np.pi
list_size = 6

SHIFT_LEFT = 0
SHIFT_RIGHT = 1
SHIFT_MIDDLE = 2
SHIFT_NONE = 3


def calc_with_lagrange(xs, ys, x):
    assert xs.size == ys.size, "Incorrect polynomial table\n"
    result = 0.0
    for i in range(xs.size):
        term = 1
        for j in range(xs.size):
            if i != j:
                term *= (x - xs[j]) / (xs[i] - xs[j])
        result += term * ys[i]
    return result


def factor(n):
    return np.prod(range(1, n + 1))


def w_fun(x, xs):
    result = 1
    for xi in xs:
        result *= x - xi
    return result


def make_list_custom(double_count=False):
    p = np.pi
    custom = [p/10, p/8, p/6, p/3, p/2, 7*p/8]
    if not double_count:
        return custom
    else:
        result = custom
        for i in range(1, len(custom)):
            result.insert(2 * i - 1, (custom[i - 1] + custom[i]) / 2)
        return result


def make_shifted_list(shift_type, double_count=False):
    left = left_bound
    right = right_bound
    if shift_type == SHIFT_LEFT:
        right /= 2
    elif shift_type == SHIFT_RIGHT:
        left = right / 2
    elif shift_type == SHIFT_MIDDLE:
        middle = (left + right) / 2
        left = (left + middle) / 2
        right = (right + middle) / 2
    elif shift_type == SHIFT_NONE:
        pass
    else:
        assert False, "incorrect shifting type"
    return np.linspace(left, right, 2 * list_size if double_count else list_size)


def process_function(fun, xs):
    format_str = "{:+15.5f}"
    n_out = 20
    ys = fun(xs)

    print "\t\tx_k = \t", "\t   f(x_k) = \t", "\tL(x_k) = \t", "\tdelta = \t", "\t   A(rate)=\t"
    for x in np.linspace(left_bound, right_bound, n_out):
        f_res = f(x)
        lagrange_res = calc_with_lagrange(xs, ys, x)

        print format_str.format(x), format_str.format(f_res), format_str.format(lagrange_res), \
            format_str.format(np.abs(lagrange_res - f_res))


def main():
    process_function(f, make_shifted_list(SHIFT_NONE))



main()

