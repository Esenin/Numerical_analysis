# coding=utf-8
from numpy.core.umath import sign

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
    while not np.all(less_xs(np.abs(np.diagonal(lam2 - lam1)))):
        lam1 = lam2
        xnew = np.dot(m, x)
        lam2 = xnew / x
        x = xnew
        k += 1
        print '\t', k, '\t', f.format(lam2[0, 0]), '\t', f.format(lam2[1, 0]), '\t', f.format(lam2[2, 0])
    return np.max(lam2)


def find_max_non_diag_hardcode(a):
    a01 = a[0, 1]
    a02 = a[0, 2]
    a12 = a[1, 2]
    maxx = np.max(np.abs(np.array([a01, a02, a12])))

    i = 0
    j = 0
    if maxx == np.abs(a01):
        j = 1
    elif maxx == np.abs(a02):
        j = 2
    else:
        i = 1
        j = 2
    return i, j


def jacob_spinner_matrix(m):
    print '\nМетод вращения Якоби:\n'
    a = m
    k = 0
    a_old = np.identity(3)

    print "A(k =", 0, ") \n"
    print_matrix(a)
    counter = 0
    while not np.all(less_xs(np.abs(np.diagonal(a - a_old)))):
        a_old = a
        k += 1
        i, j = find_max_non_diag_hardcode(a)
        maxx = a[i, j]
        phi = 0.5 * np.arctan(2 * maxx / (a[i, i] - a[j, j])) if (a[i, i] - a[j, j]) > 1e-15 else np.pi / 4

        t = np.identity(3).astype(np.float64)
        t[i, i] = np.cos(phi)
        t[j, j] = np.cos(phi)
        t[i, j] = -np.sin(phi)
        t[j, i] = np.sin(phi)

        a = np.dot(np.dot(t.T, a_old), t)

        if np.min(np.abs(a)) < 1e-13:
            counter += 1
            print "---------------------------------------------------------------------------"
            print "k = ", counter, "\t max elem = ", maxx, "\tphi = ", phi
            print "T = "
            print_matrix(t)
            print "T^T = "
            print_matrix(t.T)
            print "A( k =", counter, ") = T^T * A( k = ", counter - 1, ") * T "
            print_matrix(a)
    return np.max(a)


def jacob_spinner_equations(m):
    print '\nМетод вращения Якоби на равенствах:\n'
    counter = 0
    a = m
    a_old = np.identity(3)

    print "A(k =", 0, ") \n"
    print_matrix(a)
    while not np.all(less_xs(np.abs(np.diagonal(a - a_old)))):
        a_old = a
        counter += 1
        i, j = find_max_non_diag_hardcode(a)
        maxx = a[i, j]

        d = np.sqrt( np.power(a[i, i] - a[j, j], 2) + 4 * np.power(a[i, j], 2) )
        cosfi = np.sqrt( 0.5 * (1 + np.abs(a[i, i] - a[j, j]) / d) )
        sinfi = sign(a[i, j] * (a[i, i] - a[j, j])) * np.sqrt( 0.5 * (1 - np.abs(a[i, i] - a[j, j]) / d) )

        a = a * 0
        a[i, i] = (a_old[i, i] + a_old[j, j]) / 2 + sign(a_old[i, i] - a_old[j, j]) * d / 2
        a[j, j] = (a_old[i, i] + a_old[j, j]) / 2 - sign(a_old[i, i] - a_old[j, j]) * d / 2

        for k in range(a.shape[0]):
            for l in range(a.shape[1]):
                if (k != i) and (k != j) and (l != i) and (l != j):
                    a[k, l] = a_old[k, l]
                elif (k != i) and (k != j):
                    a[k, i] = cosfi * a_old[k, i] + sinfi * a_old[k, j]
                    a[i, k] = a[k, i]
                    a[k, j] = -sinfi * a_old[k, i] + cosfi * a_old[k, j]
                    a[j, k] = a[k, j]


        print "---------------------------------------------------------------------------"
        print "k = ", counter, "\t max elem = ", maxx
        print "A( k =", counter, ") = T^T * A( k = ", counter - 1, ") * T "
        print_matrix(a)
    return np.max(a)


def main(a_matrix):
   #a_ext = np.append(a_matrix, b_vect, axis=1)

    res1 = simple_iteration_method(a_matrix)
    res2 = jacob_spinner_equations(a_matrix)
    print "Cтепенной: ", res1, "Якоби: ",  res2
    print "delta = ", np.abs(res2 - res1)


main(A_coef)
