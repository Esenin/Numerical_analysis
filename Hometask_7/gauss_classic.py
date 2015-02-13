# coding=utf-8

__author__ = 'esenin'

'''
    Gauss classic lin.eq solver with checksum
    Ax = b
'''

import numpy as np
import pprint
import scipy
import scipy.linalg

# ######################################################################################################################
#
#   CONSTANTS
#
#######################################################################################################################
A_coef = np.matrix('3.15  1.18 2.04  0     1.32;'
                   '0     2.14 1.71  0     2.01;'
                   '-1.11 0    3.16  -3.16 1.14;'
                   '0     0    -2.73 4.20  1.50;'
                   '0     0    0     2.43 -4.51')

b_vector = np.matrix('-9.28; 1.67; -6.09; 12.66; -6.59')

format_s = "{:+8.5}"
EPS = 10 ** (-5)
#######################################################################################################################
#
#   HELPERS
#
#######################################################################################################################


def print_matrix(m):
        out = ''
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                out += '\t' + format_s.format(m[i, j])
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


def triangulate(a_m):
    rows = a_m.shape[0]
    m = 0
    for k in range(1, rows):
        for j in range(k, rows):
            m = a_m[j, k - 1] / a_m[k - 1, k - 1]
            for i in range(rows):
                a_m[j, i] = a_m[j, i] - m * a_m[k - 1, i]
    return a_m


def reversal(a_m):
    rows = a_m.shape[0]
    cols = a_m.shape[1]
    x = np.zeros((1, 5))
    for i in reversed(range(rows)):   # aka range(rows)[::-1]
        for j in range(i + 1, rows):
            x[0, i] += a_m[i, j] * x[0, j]
        x[0, i] = (a_m[i, rows] - x[0, i]) / a_m[i, i]
    return x


def main(a_matrix, b_vect):
    a_ext = np.append(a_matrix, b_vect, axis=1)
    a_ext_signed = np.append(a_ext, calc_check_vector(a_ext), axis=1)
    print '\nРасширенная матрица коэффицентов с контрольной суммой:'
    print_matrix(a_ext_signed)
    #print_matrix(triangulate(a_ext_signed))
    p, l, u = scipy.linalg.lu(a_ext_signed)
    print '\nПреобразованная матрица:'
    print_matrix(u)
    assert is_valid_checksum(u), 'Контрольные суммы не совпадают'
    print '\nРешение x* ='
    print reversal(u[:, :-1])




main(A_coef, b_vector)









