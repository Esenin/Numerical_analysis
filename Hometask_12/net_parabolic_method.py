# coding=utf-8
__author__ = 'esenin'

'''
    Метод сеток для решения уравнения параболического типа

'''

import numpy as np
# from sympy import *
# from sympy.mpmath import ln

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


f = lambda x, t: -1. * (x * (x - 1)) / (10. + t)**2 + 2 / (10.0 + t)

# u(0, t)
psi_1 = lambda t: np.exp(-t/4.0)
# u(1,t)
psi_2 = lambda t: np.exp(-t/4.) * np.cos(1./2.0)
# phi(x, 0)
phi = lambda x: np.cos(x / 2.) + x * (1 - x) / 10.0

solution = lambda x, t: np.exp(-0.25*t) * np.cos(x/2.) + x*(1-x)/10+t


def explicit(n=10):
    print "По явной схеме"
    l = 0
    r = 1
    h = (r - l) / float(n)
    m = 500
    tau = h**2 / 5
    t_max = tau * m
    print "h = ", h, "\ttau = ", tau, "\tm = ", m

    res = np.zeros((m + 2, n + 1), dtype=np.float64)
    res[0, 0] = psi_1(0)
    for i in range(1, n):
        res[0, i] = phi(l + i * h)
    res[0, n] = psi_2(0)

    for k in range(1, m + 1):
        t_k = k * tau
        for i in range(n + 1):
            x_i = l + i * h
            if x_i == l:
                res[k, i] = psi_1(t_k)
            elif x_i == r:
                res[k, i] = psi_2(t_k)
            else:
                res[k, i] = tau/h**2 * res[k-1, i+1] + (1 - 2*tau/h**2)*res[k-1, i] + tau/h**2*res[k-1, i-1] + tau*f(x_i, t_k)
    if n <= 20:
        print_matrix(res)

    return res[m, :], t_max


def main():
    #t = 500 tau, x from 0 to 1
    n = 100
    h = 1. / n
    solved, t_up = explicit(n)
    correct = [solution(0 + i * h, t_up) for i in range(n)]

    for i in range(n):
        print "delta = ", (solved[i] - correct[i])




main()