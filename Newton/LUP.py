import numpy

def LUP(A, P, Q, n, m):
    C = [[0] * m for i in range(n)]
    P=numpy.identity(10)
    Q=numpy.identity(10)
    for i in range(n):
        for j in range(m):
            C[i][j] = A[i][j]

    for i in range(n):
        pivotValue = 0
        pivot1 = -1
        pivot2 = -1

        for j in range(i, n):
            for k in range(i, m):
                if abs(C[j][k]) > pivotValue:
                    pivotValue = abs(C[j][k])
                    pivot1 = j
                    pivot2 = k

        if pivotValue != 0:
            for k in range(m):
                P[pivot1][k], P[i][k] = P[i][k], P[pivot1][k]
                C[pivot1][k], C[i][k] = C[i][k], C[pivot1][k]
            for j in range(n):
                Q[j][pivot2], Q[j][i] = Q[j][i], Q[j][pivot2]
                C[j][pivot2], C[j][i] = C[j][i], C[j][pivot2]
            for j in range(i + 1, n):
                C[j][i] /= C[i][i]
                for k in range(i + 1, m):
                    C[j][k] -= C[j][i] * C[i][k]

    return C,P,Q


def U(G, n, m):
    U = [[0] * m for i in range(n)]
    for i in range(n):
        for j in range(m):
            if j < i:
                U[i][j] = 0
            else:
                U[i][j] = G[i][j]

    return U


def L(G, n, m):
    L = [[0] * n for i in range(n)]
    for i in range(n):
        for j in range(n):
            if j > i:
                L[i][j] = 0
            elif j == i:
                L[i][j] = 1
            else:
                L[i][j] = G[i][j]
    return L


def Podstanovka(L, U, n, b, m):
    rang = n
    x = [[0] for i in range(n)]
    y = [[0] for i in range(n)]

    for i in range(rang):
        s = sum(L[i][j] * y[j][0] for j in range(i))
        y[i][0] = b[i][0] - s

    for i in reversed(range(rang)):
        s = sum(U[i][j] * x[j][0] for j in range(i + 1, rang))
        x[i][0] = (y[i][0] - s) / U[i][i]

    return x
