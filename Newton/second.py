import math
import numpy as np
import LUP
from math import *
import time
import Jacob

print('Метод решения:\n1. Полный\n2. Модифицированный\n3. Гибридный\n ')
option = int(input())

if option == 2:
    print('Введите номер итерации\n')
    nomber = int(input())
elif option == 3:
    print('Введите количество итераций\n')
    nomber = int(input())

start_time = time.time()

X = [[0.5],
     [0.5],
     [1.5],
     [-1],
     [-0.2],
     [1.5],
     [0.5],
     [-0.5],
     [1.5],
     [-1.5]]

X0 = [[1] for i in range(10)]

T = [[0] * 10 for i in range(10)]

iteration = 0
operation = 0
L = np.identity(10)
U = np.identity(10)
P = np.identity(10)
Q = np.identity(10)

while (max((X[i][0] - X0[i][0]) ** 2 for i in range(10)) > 10e-10):

    X0 = np.copy(X)


    F = Jacob.func(X)
    operation += 10

    F0 = [[-F[0][0]],
          [-F[1][0]],
          [-F[2][0]],
          [-F[3][0]],
          [-F[4][0]],
          [-F[5][0]],
          [-F[6][0]],
          [-F[7][0]],
          [-F[8][0]],
          [-F[9][0]]]

    if option == 2:
        if iteration < nomber:
            T = Jacob.jacob(X, T)
            operation += 67
            G,P,Q = LUP.LUP(T, P, Q, 10, 10)
            U = LUP.U(G, 10, 10)
            L = LUP.L(G, 10, 10)
            operation += 667
    elif option == 3:
        if iteration % nomber == 0 or iteration>=4:
            T = Jacob.jacob(X, T)
            operation += 67
            G,P,Q = LUP.LUP(T, P, Q, 10, 10)
            U = LUP.U(G, 10, 10)
            L = LUP.L(G, 10, 10)
            operation += 667
    else:
        T = Jacob.jacob(X, T)
        operation += 67
        G,P,Q = LUP.LUP(T, P, Q, 10, 10)
        U = LUP.U(G, 10, 10)
        L = LUP.L(G, 10, 10)
        operation += 667


    S = np.dot(P, F0)
    S = LUP.Podstanovka(L, U, 10, S, 10)
    operation += 100
    S = np.dot(Q, S)
    X = X0 + S
    print(X)
    print()
    iteration += 1

print(iteration, ' iterations')
print(operation, ' operations')
print("--- %s seconds ---" % (time.time() - start_time))

F1 = Jacob.func(X)

print(F1)
