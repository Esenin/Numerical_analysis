# coding=utf-8

__author__ = 'esenin'

'''
    Roots method for lin.eq with checksum
    Ax = b, where A = S*S'

                   ( a[0, 0]               | f[0]   | f~[0]   )
    A_ext_signed = (        ...            | ...    | ...     )
                   (           a[n-1, n-1] | f[n-1] | f~[n-1] )

                   (s[0, 0]                | k[0]   | k~[0]   )
    S_ext_signed = (       ...             | ...    | ...     )
                   (           s[n-1, n-1] | k[n-1] | k~[n-1] )
'''

import numpy as np
import scipy
import scipy.linalg

# ######################################################################################################################
#
#   CONSTANTS
#
#######################################################################################################################
A_coef = np.matrix('-3.045  2.113  2.113  -2.113 ;'
                   '2.113   3.045  1.798   1.798 ;'
                   '2.113   1.798  2.563   1.798 ;'
                   '-2.113  1.798  1.798  -3.074')


b_vector = np.matrix('4.475; 13.361; 10.385; 6.040')

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


def calc_check_vector(a_m):
    return a_m.sum(axis=1)


def is_valid_checksum(a_vo):
    print 'Пересчет вектора контрольной суммы:'
    checksum = np.matrix(calc_check_vector(a_vo[:, :-1])).T
    print checksum
    embedded_sum = a_vo[:, -1:]
    assert checksum.shape == embedded_sum.shape, 'Не совпадает размерность контрольных векторов'
    return all([abs(checksum[i, 0] - embedded_sum[i, 0]) < EPS for i in range(checksum.shape[0])])


def choleskiy_decomposition(a_m):
    a = a_m.astype(np.complex128)
    rows = a_m.shape[0]
    s_m = a.copy()
    s_m.fill(0+0.j)
    for i in range(rows):
        for j in range(i):
            temp = 0.0
            for k in range(j):
                temp += s_m[i, k] * s_m[j, k]

            s_m[i, j] = (a[i, j] - temp) / s_m[j, j]

        temp2 = a[i, i]
        for k in range(i):
            temp2 -= s_m[i, k] * s_m[i, k]
        s_m[i, i] = np.sqrt(temp2)
    return s_m


def k_recalc(a_em, s):
    rows = a_em.shape[0]
    k = np.zeros((4, 1)).astype(np.complex128)
    k[0, 0] = a_em[0, 4] / s[0, 0]
    for i in range(1, rows):
        temp = a_em[i, 4]
        for l in range(i):
            temp -= s[l, i] * k[l, 0]
        k[i, 0] = temp / s[i, i]
    return k


def reversal(s, k):
    rows = s.shape[0]
    x = np.zeros((4, 1)).astype(np.complex128)

    x[3, 0] = k[3, 0] / s[rows - 1, rows - 1]

    for i in reversed(range(rows)):
        temp = k[i, 0]
        for l in range(i + 1, rows):
            temp -= s[i, l] * x[l, 0]
        x[i, 0] = temp / s[i, i]
    return x


def main(a_matrix, b_vect):
    a_ext = np.append(a_matrix, b_vect, axis=1)
    a_ext_signed = np.append(a_ext, calc_check_vector(a_ext), axis=1)
    print '\nРасширенная матрица коэффицентов с контрольной суммой:'
    print_matrix(a_ext_signed)

    s_m = choleskiy_decomposition(A_coef)
    k = k_recalc(a_ext.astype(np.complex128), s_m)
    s_ext = np.append(s_m, k, axis=1)
    s_ext_signed = np.append(s_ext, calc_check_vector(s_ext), axis=1)
    print '\nПреобразованная матрица:'
    print_matrix(s_ext_signed, is_complex=True)
    #assert is_valid_checksum(s_m), 'Контрольные суммы не совпадают'
    print '\nРешение x* ='
    print reversal(s_m, k)




main(A_coef, b_vector)