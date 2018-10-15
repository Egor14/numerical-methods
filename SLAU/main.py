import math
import numpy as np
from Matrix import Matrix
import random
import copy

Mat = Matrix()
start=1
while start==1:
    print('Введите размерность')
    n = int(input())
    m = int(input())

    print('1. Ввести матрицу\n2. Сгенерировать матрицу')
    option = int(input())
    if option == 1:
        A = [[float(j) for j in input().split()] for i in range(n)]
        print('Введите вектор b')
        b = []
        b0 = []
        for i in range(n):
            x = float(input())
            b.append(x)
            b0.append(x)
    else:
        print('1. Любая матрица\n2. Матрица с диагональным преобладанием\n3. Положительно определенная матрица')
        number=int(input())
        A = np.random.randint(-100, 100, (n, m))
        if number ==2:
            h = min(n, m)
            for i in range(h):
                A[i][i] = sum(math.fabs(A[i][j]) for j in range(m))
        elif number==3:
            A=np.dot(A,np.transpose(A))
            for i in range(n):
                for j in range(m):
                    A[i][j]=A[i][j]/100
        Mat.MatPrint(A, n, m)
        b = []
        b0 = []
        for i in range(n):
            x = random.randint(-100, 100)
            b.append(x)
            b0.append(x)
        print(b)

    print('\n1. LUP\n2. QR\n3. Якоби\n')
    option = int(input())
 
    if option == 1:
        P = np.identity(n)
        Q = np.identity(m)

        G = Mat.LUP(A, P, Q, n, m)

        U = Mat.U(G, n, m)
        L = Mat.L(G, n, m)

        print()
        print('C=')
        Mat.MatPrint(G, n, m)

        # Вывод матриц L, U и проверка
        print()
        print('U=')
        Mat.MatPrint(U, n, m)
        print()
        print('L=')
        Mat.MatPrint(L, n, n)
        # M1 = Mat.MULT(Mat.MULT(P, A, n, m),Q,n,m)
        M1 = Mat.MULT(P, A, n, m)
        M1 = Mat.MULT(M1, Q, n, m)
        M2 = Mat.MULT(L, U, n, m)
        print()
        print('P * A * Q =')
        Mat.MatPrint(M1, n, m)
        print()
        print('L * U =')
        Mat.MatPrint(M2, n, m)
        print()
        print('Q')
        Mat.MatPrint(Q, n, n)
        print('P')
        Mat.MatPrint(P, n, n)

        if n == m:
            # Определитель и ранг
            det = Mat.Det(U, n)
            print('Det = ', '%.3f' % det)
            print()
            rang = Mat.Rang(U, n, m)
            print('Rang = ', rang)
            print()

            # Обратная матрица с проверкой
            if det != 0:
                print('Reversed matrix =')
                R = Mat.RevMatrix(A, L, U, n, m)
                R=np.dot(Q,R)
                R=np.dot(R,P)
                R = np.transpose(R)
                Mat.MatPrint(R, n, n)
                print()
                print('A * R = ')
                print(Mat.MatPrint(Mat.MULT(A, R, n, n), n, n))
                print()

        # Решение системыс проверкой
        b = Mat.VecMult(P, b, n, n)
        b = Mat.Podstanovka2(L, U, n, b, m, 1, Q)
        b = Mat.VecMult(Q, b, m, m)
        print('Решение', b)
        print()
        print('Check: A*x - b =')
        b = ((Mat.VecMult(A, b, n, m)))
        for i in range(n):
            print('%.3f' % b[i])
        print('-')
        print(b0)
        print()

        # Число обусловленности матрицы
        if n == m:
            print('Число обусловленности матрицы А = ', Mat.Norma(A, n) * Mat.Norma(Mat.RevMatrix(A, L, U, n, m), n))
    elif option == 2:
        print('\nQR\n')
        Q = Mat.QR(A, n)
        R = Mat.MULT(Q, A, n, n)
        print('\nQ\n')
        Mat.MatPrint(Q, n, n)
        print('\nR\n')
        Mat.MatPrint(R, n, n)
        b = Mat.VecMult(Q, b, n, n)
        b = Mat.Podstanovka(R, R, n, b, n, 0)
        print('\nРешение :')
        for i in range(n):
            print('%.3f' % b[i],end=' ')
        print()
        print()
        b = ((Mat.VecMult(A, b, n, m)))
        for i in range(n):
            print('%.3f' % b[i])
        print('-')
        print(b0)
        print()

    else:
        C = [0] * n
        B = Mat.Jacobi(A, C, b, n)
        Mat.MatPrint(B, n, n)
        print('\n', C)
        b = Mat.First(B, C, n)
        print('\nРешение :')
        for i in range(n):
            print('%.3f' % b[i],end=' ')
        print()
        print()
        b = ((Mat.VecMult(A, b, n, m)))
        for i in range(n):
            print('%.3f' % b[i])
        print('-')
        print(b0)
        print()
        b = Mat.Seidel(B, C, n)
        print('\nРешение :')
        for i in range(n):
            print('%.3f' % b[i], end=' ')
        print()
        print()
        b = ((Mat.VecMult(A, b, n, m)))
        for i in range(n):
            print('%.3f' % b[i])
        print('-')
        print(b0)
        print()
    start=int(input())


