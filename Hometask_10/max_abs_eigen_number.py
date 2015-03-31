# coding=utf-8

__author__ = 'esenin'

'''
    Степенной метод и метод Якоби поиска макс. по модулю собственные числа и вектор

'''


import numpy as np
#import scipy
#import scipy.linalg

# ######################################################################################################################
#
#   CONSTANTS
#
#######################################################################################################################
A_coef = np.matrix('.72     .45     .38 ;'
                   '.45     .91     .56 ;'
                   '.38     .56     .12 ')


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


def less_eps(x):
    return x < EPS

less_xs = np.vectorize(less_eps)


def simple_iteration_method(m):
    print '\nСтепенной метод:\n'
    print '\tk\t\tx1\t\t\tx2\t\t\tx3'

    f = format_s
    x = np.ones((3, 1)).astype(np.float64)
    lam1 = np.ones((3, 1)).astype(np.float64)
    lam2 = lam1 * 2
    k = 0
    while not np.any(less_xs(np.abs(np.diagonal(lam2 - lam1)))):
        lam1 = lam2
        xnew = np.dot(m, x)
        lam2 = xnew / x
        x = xnew
        k += 1
        print '\t', k, '\t', f.format(lam2[0, 0]), '\t', f.format(lam2[1, 0]), '\t', f.format(lam2[2, 0])


def jacob_spinner(m):
    print '\nМетод вращения Якоби:\n'
    a = m
    k = 0
    a_old = np.identity(3)
    while not np.all(less_xs(np.abs(a - a_old))):
        k += 1
        a01 = a[0, 1]
        a02 = a[0, 2]
        a12 = a[1, 2]
        maxx = np.max(np.array([a01, a02, a12]))
        if np.abs(maxx) < EPS:
            break
        i = 0
        j = 0
        if maxx == a01:
            j = 1
        elif maxx == a02:
            j = 2
        else:
            i = 1
            j = 2
        phi = 0.5 * np.arctan(2 * maxx / (a[i, i] - a[j, j])) if (a[i, i] - a[j, j]) > EPS else np.pi / 4

        T = np.identity(3).astype(np.float64)
        T[i, i] = np.cos(phi)
        T[j, j] = T[i, i]
        T[j, i] = np.sin(phi)
        T[i, j] = -T[j, i]

        a = np.dot(T.T, np.dot(a, T))
        print "\nk = ", k, "\nT = \n"
        print_matrix(T)
        print "A(k) = \n"
        print_matrix(a)


def main(a_matrix):
   #a_ext = np.append(a_matrix, b_vect, axis=1)

    simple_iteration_method(a_matrix)
    jacob_spinner(a_matrix)


main(A_coef)
