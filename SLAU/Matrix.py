import math
import numpy as np

class Matrix:
    def LUP(self, A, P, Q, n, m):
        C = [[0] * m for i in range(n)]
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


        return C

    def U(self, G, n, m):
        U = [[0] * m for i in range(n)]
        for i in range(n):
            for j in range(m):
                if j < i:
                    U[i][j] = 0
                else:
                    U[i][j] = G[i][j]
        return U

    def L(self, G, n, m):
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



    def MULT(self, A, B, n, m):
        M = [[0] * m for i in range(n)]
        for i in range(n):
            for j in range(m):
                M[i][j] = 0
                for k in range(n):
                    M[i][j] += A[i][k] * B[k][j]
                if abs(M[i][j])<10**(-10):
                    M[i][j]=0
        return M

    def VecMult(self, P,b,n,m):
        x=[0]*n
        for i in range(n):
            for j in range(m):
                x[i]+=P[i][j]*b[j]
            if abs(x[i])<10**(-10):
                x[i]=0
        return x


    def MatPrint(self, A, n, m):
        for i in range(n):
            print()
            for j in range(m):
                print('%.3f' % A[i][j], end=' ')
        print()



    def Det(self, U, n):
        det1=1
        for i in range(n):
            det1*=U[i][i]
        return det1

    def Rang(self,U,n,m):
        # rang=n
        # while math.fabs(U[rang-1][m-1])<(10e-14):
        #     rang-=1
        rang = 0
        i=0
        j=0
        while i<n:
            j=i
            while j<m:
                if U[i][j]!=0:
                    j=m+1
                j+=1
            if j!=m:
                rang+=1
            i+=1


        return rang


    def Podstanovka(self, L,U,n,b,m,k):
        rang=self.Rang(U,n,m)
        y=b
        x=m*[0]

        if k==1:
            y = m * [0]
            for i in range(rang):
                s = sum(L[i][j] * y[j] for j in range(i))
                y[i] = b[i] - s



        for i in reversed(range(rang)):
            s=sum(U[i][j]*x[j] for j in range(i+1, rang))
            x[i]=(y[i]-s)/U[i][i]


        return x

    def Podstanovka2(self, L,U,n,b,m,k,Q):
        rang=self.Rang(U,n,m)
        y=b
        x=m*[0]

        if k==1:
            y = m * [0]
            for i in range(rang):
                s = sum(L[i][j] * y[j] for j in range(i))
                y[i] = b[i] - s


        for i in reversed(range(rang)):
            s=sum(U[i][j]*x[j] for j in range(i+1, rang))
            x[i]=(y[i]-s)/U[i][i]


        return x

    def RevMatrix(self,A,L,U,n,m):
        R=[]
        E=[0]*n
        for i in range(n):
            if i==0:
                E[i]=1
            else:
                E[i-1]=0
                E[i]=1

            R.append( self.Podstanovka(L,U,n,E,m,1))
        return R

    def Norma(self,A,n):
        norm=0
        for i in range(n):
            for j in range(n):
                norm+=(A[i][j])*A[i][j]
        return math.sqrt(norm)

    def VecRaz(self, x ,y, n):
        #m=abs(max(x[i]-y[i] for i in range(n)))
        m=math.sqrt(sum((x[i]-y[i])**2 for i in range(n)))
        return m

    def Jacobi(self, A, C, b, n):
        B = [[0] * n for i in range(n)]
        for i in range(n):
            C[i]=b[i]/A[i][i]
            for j in range(n):
                if i != j:
                    B[i][j]=-A[i][j]/A[i][i]
        return B

    def Seidel(self,B,C,n):
        x0=[0]*n
        x1=[1]*n
        k=0
        while self.VecRaz(x0,x1,n)>10**(-13):
            x1=np.copy(x0)
            k+=1
            for i in range(n):
                x0[i]=sum(B[i][j]*x0[j] for j in range(n))
                x0[i]+=C[i]

        print('\nМетод Зейделя получил решение за ', k, ' итераций')

        return x0

    def First(self,B,C,n):
        x0=[0]*n
        x1=[1]*n
        x2=[0]*n
        k=0
        while self.VecRaz(x0,x1,n)>10**(-13):
            x1=np.copy(x0)
            k+=1
            for i in range(n):
                x0[i]=sum(B[i][j]*x2[j] for j in range(n))
                x0[i]+=C[i]
            x2=np.copy(x0)

        print('\nПервый метод получил решение за ',k,' итераций')

        return x0


    def Raz(self, A, B, n):
        C = [[0] * n for i in range(n)]
        for i in range(n):
            for j in range(n):
                C[i][j]=(A[i][j]-B[i][j])

        return C

    def QR(self, A, n):
        w0 = [[0] for i in range(n)]
        w1 = [[0] * n for i in range(1)]
        E=np.identity(n)
        Q=np.identity(n)
        j=0
        for i in range(n-1):
            for t in range(i, n):
                w0[t][0] = A[t][i]
                w1[0][t] = A[t][i]
            alpfa=math.sqrt(sum(A[k][i]**2 for k in range(j, n)))
            al=sum(A[k][i]**2 for k in range(j+1, n))
            al=(al+(A[j][i]-alpfa)**2)
            w0[i][0]-=alpfa
            w1[0][i]-=alpfa
            if al==0:
                al=1
            V0=self.Raz(E, self.Mult(w0,w1,n, al),n)
            print('\nV0\n')
            V0=self.V(V0,n,i)
            self.MatPrint(V0,n,n)
            Q=self.MULT(V0,Q,n,n)
            A=self.MULT(V0,A,n,n)
            j+=1
            w0[i][0] = 0
            w1[0][i] = 0

        return Q


    def V(self,V0,n,i):
        if i>0:
            for j in range(n):
                for k in range(n):
                    if k<i:
                        V0[k][j]=0
                        V0[j][k]=0
                        V0[k][k]=1
            V0[i - 1][i - 1] = 1

        return V0


    def Mult(self, w0, w1, n, al):
        H = [[0] * n for i in range(n)]
        for i in range(n):
            for j in range(n):
                H[i][j]=(w0[i][0]*w1[0][j])*2/al

        return H