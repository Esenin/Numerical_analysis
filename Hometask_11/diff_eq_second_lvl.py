# coding=utf-8
__author__ = 'esenin'

'''
    Метод прогонки

'''

import numpy as np
from sympy import *
#from sympy.mpmath import ln

# ######################################################################################################################
#
#   CONSTANTS
#
#######################################################################################################################

format_s = "{:+8.5}"
format_cs = "{:+15.5}"
EPS = 10 ** (-5)
#######################################################################################################################
#
#   HELPERS
#
#######################################################################################################################


def print_matrix(m, is_complex=False):
    f = format_cs if is_complex else format_s
    out = ''
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            out += '\t' + f.format(m[i, j])
        out += '\n'
    print out

#######################################################################################################################
#
#   GENERAL
#
#######################################################################################################################
x = Symbol('x')
model_y = x**3 + 2*x - 1


def thomas_algo(alpha, betta, fx, imodel = False, n=100):
    x0 = 0
    xn = 1
    h = float(xn - x0) / n

    xs = [x0 + i * h for i in range(0, n + 1)]
    ysym = [Symbol('y' + str(i)) for i in range(0, n + 1)]

    px = (2 + ln(x + 1))
    qx = -(1 + x**2)
    # fx = 1 - x

    a = 1 + px * h/2
    b = 2 - qx * h**2
    c = 1 - px * h/2
    g = fx * h**2

    kappa1 = ((b - 4*a) / (c - 2 * a * h * alpha - 3 * a)).subs(x, xs[1])
    nu1    = (g / (c - 2 * a * h * alpha - 3 * a)).subs(x, xs[1])

    kappa2 = ((4*c - b) / (3 * c - 2 * c * h * betta - a)).subs(x, xs[n-1])
    nu2    = (-g / (3 * c - 2 * c * h * betta - a)).subs(x, xs[n-1])

    us = [0] * (n)
    vs = [0] * (n)
    us[0] = kappa1
    vs[0] = nu1
    print "\t\t y[i] = u[i] y[i+1] + v[i]"
    for i in range(1, n):
        us[i] = (a / (b - c * us[i-1])).subs(x, xs[i])
        vs[i] = ((c * vs[i-1] - g) / (b - c * us[i-1])).subs(x, xs[i])
        print "y[", i, "] = {:+8.3}".format(us[i]), "\ty[", i+1, "] + {:+8.3}".format(vs[i])


    ys = [0] * (n+1)
    ys[n] = float((nu2 + kappa2 * vs[n-1]) / (1 - kappa2 * us[n-1]))

    for i in range(n-1, -1, -1):
        ys[i] = solve(ysym[i] - us[i] * ys[i+1] - vs[i], ysym[i])[0]

    print '\t x\t\t\ty', "model y(x_i) =" if imodel else ""
    for (xi, yi) in zip(xs, ys):
        print format_s.format(xi), "\t", format_s.format(yi), "\t ",  model_y.subs(x, xi) if imodel else ""


def main():
    print "Модельное решение:"
    ffun = 5 + 4 * x + 7 * x**2 - 3 * x**3 - x**5 + (2 + 3 * x**2) * ln(1 + x)
    thomas_algo(-2., 5/2., ffun, True)

    print "-" * 110
    print "Основная задача:"
    ffun_m = 1 - x
    thomas_algo(1.2, -0.4, ffun_m, False)




main()
