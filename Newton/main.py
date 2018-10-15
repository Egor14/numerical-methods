import math
import numpy as np
import LUP
from math import *
import time

X = [[0],
     [1]]

X0 = [[2],
      [1]]

while (max((X[i][0] - X0[i][0]) ** 2 for i in range(2)) > 10e-4):
    X0 = np.copy(X)
    P = np.identity(2)
    Q = np.identity(2)

    F = [[math.sin(X[0][0]) + 2 * X[1][0] - 2],
         [X[0][0] + math.cos(X[1][0] - 1) - 0.7]]

    F0 = [[-F[0][0]],
          [-F[1][0]]]

    T = [[math.cos(X[0][0]), 2],
         [1, -math.sin(X[1][0] - 1)]]

    G = LUP.LUP(T, P, Q, 2, 2)
    U = LUP.U(G, 2, 2)
    L = LUP.L(G, 2, 2)
    S = np.dot(P, F0)
    S = LUP.Podstanovka(L, U, 2, S, 2)
    S = np.dot(Q, S)
    X = X0 + S
    print(X)
    print()

