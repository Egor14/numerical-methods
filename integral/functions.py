import numpy as np
from math import *
import scipy.integrate as integrate


def f(x):
    return cos(1.5 * x) * exp(2 * x / 3) + 3 * sin(5.5 * x) * exp(-2 * x) + 2


def m0(a, b):
    #return integrate.quad(lambda x: 1 / (x - 2.5) ** (2 / 7), a, b)[0]
     return 1.4 * (b - 2.5) ** (5 / 7) - 1.4 * (a - 2.5) ** (5 / 7)


def m1(a, b):
    #return integrate.quad(lambda x: x / (x - 2.5) ** (2 / 7), a, b)[0]
     return 7 / 24 * (b - 2.5) ** (5 / 7) * (2 * b + 7) - 7 / 24 * (a - 2.5) ** (5 / 7) * (2 * a + 7)


def m2(a, b):
   # return integrate.quad(lambda x: x ** 2 / (x - 2.5) ** (2 / 7), a, b)[0]
     return 7 / 456 * (b - 2.5) ** (5 / 7) * (24 * b ** 2 + 70 * b + 245) - 7 / 456 * (a - 2.5) ** (5 / 7) * (
        24 * a ** 2 + 70 * a + 245)


def m3(a, b):
    #return integrate.quad(lambda x: x ** 3 / (x - 2.5) ** (2 / 7), a, b)[0]
    return 7 / 7904 * (b - 2.5) ** (5 / 7) * (304 * b ** 3 + 840 * b ** 2 + 2450 * b + 8575) - 7 / 7904 * (
                                                                                                               a - 2.5) ** (
                                                                                                               5 / 7) * (
                                                                                                    304 * a ** 3 + 840 * a ** 2 + 2450 * a + 8575)


def m4(a, b):
    #return integrate.quad(lambda x: x ** 4 / (x - 2.5) ** (2 / 7), a, b)[0]
     return 7 / 130416 * (b - 2.5) ** (5 / 7) * (
         3952 * b ** 4 + 10640 * b ** 3 + 29400 * b ** 2 + 85750 * b + 300125) - 7 / 130416 * (a - 2.5) ** (5 / 7) * (
         3952 * a ** 4 + 10640 * a ** 3 + 29400 * a ** 2 + 85750 * a + 300125)


def m5(a, b):
    #return integrate.quad(lambda x: x ** 5 / (x - 2.5) ** (2 / 7), a, b)[0]
     return 7 / 10433280 * (b - 2.5) ** (5 / 7) * (
         260832 * b ** 5 + 691600 * b ** 4 + 1862000 * b ** 3 + 5145000 * b ** 2 + 15006250 * b + 52521875) - 7 / 10433280 * (
                                                                                                                                 a - 2.5) ** (
                                                                                                                                 5 / 7) * (
                                                                                                                  260832 * a ** 5 + 691600 * a ** 4 + 1862000 * a ** 3 + 5145000 * a ** 2 + 15006250 * a + 52521875)


def IKF(a, b):
    x = [a, ((2*a + b) / 3), b]
    u = [0] * 3
    u[0] = m0(a, b)
    u[1] = m1(a, b)
    u[2] = m2(a, b)
    X = [[1] * 3 for i in range(3)]
    X[1][0] = x[0]
    X[1][1] = x[1]
    X[1][2] = x[2]
    X[2][0] = x[0] ** 2
    X[2][1] = x[1] ** 2
    X[2][2] = x[2] ** 2
    A = np.linalg.solve(X, u)
    S = 0
    for i in range(3):
        S += A[i] * f(x[i])
    return S


def Gauss(a, b):
    u = [0] * 6
    u[0] = m0(a, b)
    u[1] = m1(a, b)
    u[2] = m2(a, b)
    u[3] = m3(a, b)
    u[4] = m4(a, b)
    u[5] = m5(a, b)
    A = [[0] * 3 for i in range(3)]
    for i in range(3):
        for j in range(3):
            A[i][j] = u[i + j]
    y = [0] * 3
    y[0] = -u[3]
    y[1] = -u[4]
    y[2] = -u[5]
    x = np.linalg.solve(A, y)
    x = polynom(x)
    for i in range(3):
        if x[i] > 4.3 or x[i] < 2.5:
            print("error")
    X = [[1] * 3 for i in range(3)]
    X[1][0] = x[0]
    X[1][1] = x[1]
    X[1][2] = x[2]
    X[2][0] = x[0] ** 2
    X[2][1] = x[1] ** 2
    X[2][2] = x[2] ** 2
    A = np.linalg.solve(X, [u[0], u[1], u[2]])
    S = 0
    for i in range(3):
        S += A[i] * f(x[i])
    return S


def runge(S1, S2, L, m):
    return abs((S2 - S1) / (1 - L ** (-m)))


def eitken(S1, S2, S3, L):
    return -np.log(abs((S3 - S2) / (S2 - S1))) / np.log(L)


def polynom(a):
    p = a[1] - (a[2] ** 2) / 3
    q = a[0] + 2 * (a[2] ** 3) / 27 - a[2] * a[1] / 3
    d = q ** 2 / 4 + p ** 3 / 27
    g = np.arccos((-q / 2) * (-3 / p) ** (3 / 2))
    x = [0] * 3
    x[0] = 2 * np.sqrt(-p / 3) * np.cos(g / 3) - a[2] / 3
    x[1] = 2 * np.sqrt(-p / 3) * np.cos(g / 3 + 2 * np.pi / 3) - a[2] / 3
    x[2] = 2 * np.sqrt(-p / 3) * np.cos(g / 3 - 2 * np.pi / 3) - a[2] / 3
    return x
