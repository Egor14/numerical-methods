import numpy as np
import scipy.integrate as integrate

def solveF(x):
    return 1.3*np.cos(3.5*x)*np.exp(2*x/3) + 6*np.sin(4.5*x)*np.exp((-x)/8) + 5*x
def solveP(x):
    return 1/((3.2-x)**0.25)
def solvePW(x):
    return abs((x-0.7)*(x-(0.7+3.2)/2)*(x-3.2)/((3.2-x)**0.25))
def solvePWGaus(x):
    return abs((x-0.7)**2*(x-(0.7+3.2)/2)**2*(x-3.2)**2/((3.2-x)**0.25))


def methodicLatency():
    return abs((100)* integrate.quad(lambda x: solvePW(x),0.7,3.2)[0])
def methodicLatencyGaus():
    return abs((100)* integrate.quad(lambda x: solvePWGaus(x),0.7,3.2)[0])
def analiticLatency(S):
    return abs(integrate.quad(lambda x: solveP(x)*solveF(x),0.7,3.2)[0] - S)


def solveP0(x2,x1):
    return integrate.quad(lambda x: solvePW(x),x2,x1)[0]
    #return ((4*x1)/3 - 64/15)/(16/5 - x1)**(1/4)-((4*x2)/3 - 64/15)/(16/5 - x2)**(1/4)
def solveP1(x2,x1):
    return integrate.quad(lambda x: x*solvePW(x),x2,x1)[0]
    #return ((4*x1**2)/7 + (64*x1)/105 - 4096/525)/(16/5 - x1)**(1/4)-((4*x2**2)/7 + (64*x2)/105 - 4096/525)/(16/5 - x2)**(1/4)
def solveP2(x2,x1):
    return integrate.quad(lambda x: x**2*solvePW(x),x2,x1)[0]
    #return -4/11*(16/5 - x1)**(3/4)*(x1**2 + (128*x1)/35 + 8192/525)+4/11*(16/5 - x2)**(3/4)*(x2**2 + (128*x2)/35 + 8192/525)
def solveP3(x2,x1):
    return integrate.quad(lambda x: x**3*solvePW(x),x2,x1)[0]
    #return -4/15*(16/5 - x1)**(3/4)*(x1**3 + (192*x1**2)/55 + (24576*x1)/1925 + 524288/9625)+4/15*(16/5 - x2)**(3/4)*(x2**3 + (192*x2**2)/55 + (24576*x2)/1925 + 524288/9625)
def solveP4(x2,x1):
    return integrate.quad(lambda x: x**4*solvePW(x),x2,x1)[0]
    #return  ((4*x1**5)/19 + (64*x1**4)/1425 + (16384*x1**3)/78375 + 1.14677*x1**2 + 9.78575**x1 - 125.258)/(16/5 - x1)**(1/4)-((4*x2**5)/19 + (64*x2**4)/1425 + (16384*x2**3)/78375 + 1.14677*x2**2 + 9.78575**x2 - 125.258)/(16/5 - x2)**(1/4)
def solveP5(x2,x1):
    return integrate.quad(lambda x: x**5*solvePW(x),x2,x1)[0]
    #return ((4*x1**6)/23 + (64*x1**5)/2185 + (4096*x1**4)/32775 + 0.581694*x1**3 + 3.19101*x1**2 + 27.2299*x1 - 348.543)/(16/5 - x1)**(1/4)-((4*x2**6)/23 + (64*x2**5)/2185 + (4096*x2**4)/32775 + 0.581694*x2**3 + 3.19101*x2**2 + 27.2299*x2 - 348.543)/(16/5 - x2)**(1/4)
#исправить погрешность для моментов
def rungeRule(S1,S2,L,m = 4):
    return abs((S2-S1)/(1-L**(-m)))
def eitkenRule(S1,S2,S3,L):
    return -np.log(abs((S3-S2)/(S2-S1)))/np.log(L)


def buildIKF(a,b):
    x = [a ,(a+b)/2 , b]
    u = [0]*3
    u[0] = solveP0(a,b)
    u[1] = solveP1(a,b)
    u[2] = solveP2(a,b)
    X = [[1]*3 for i in range(3)]
    X[1][0] = x[0]
    X[1][1] = x[1]
    X[1][2] = x[2]
    X[2][0] = x[0] ** 2
    X[2][1] = x[1] ** 2
    X[2][2] = x[2] ** 2
    A = np.linalg.solve(X,u)
    S = 0
    for i in range (3):
        S+= A[i]*solveF(x[i])
    return S
def buildGaus(a,b):
    u = [0]*6
    u[0] = solveP0(a,b)
    u[1] = solveP1(a,b)
    u[2] = solveP2(a,b)
    u[3] = solveP3(a,b)
    u[4] = solveP4(a,b)
    u[5] = solveP5(a,b)
    A = [[0] * 3 for i in range(3)]
    for i in range(3):
        for j in range(3):
            A[i][j]=u[i+j]
    b1 = [0]*3
    b1[0] = -u[3]
    b1[1] = -u[4]
    b1[2] = -u[5]
    x = solvePol(np.linalg.solve(A,b1))
    X = [[1]*3 for i in range(3)]
    for i in range(3):
        if x[i]>3.2 or x[i]<0.7:
            print("smthwrong")
    X[1][0] = x[0]
    X[1][1] = x[1]
    X[1][2] = x[2]
    X[2][0] = x[0] ** 2
    X[2][1] = x[1] ** 2
    X[2][2] = x[2] ** 2
    A = np.linalg.solve(X,[u[0],u[1],u[2]])
    S = 0
    for i in range (3):
        S+= A[i]*solveF(x[i])
    return S



def solvePol(a):
    p = a[1] - (a[2]**2)/3 ## -0.983152419
    q = a[0] + 2*(a[2]**3)/27 - a[2]*a[1]/3  ##0.0217037
    d = q**2/4+p**3/27
    g = np.arccos((-q/2)*(-3/p)**(3/2))
    x = [0]*3
    x[0] = 2 * np.sqrt(-p / 3) * np.cos(g / 3) - a[2]/3
    x[1] = 2 * np.sqrt(-p / 3) * np.cos(g / 3 + 2*np.pi/3) - a[2]/3
    x[2] = 2 * np.sqrt(-p / 3) * np.cos(g / 3 - 2*np.pi/3) - a[2]/3
    return x

S1 = buildIKF(0.7,3.2)
print(S1)
print(methodicLatency())
print(analiticLatency(S1))
h=2.5
maxE = 100
k = 1
maxS = 1
count = 1
while(abs(maxE)>10**-6):
    k+=1
    S = 0
    h = h/2
    a = 0.7
    b = 3.2
    count = count * 2
    z=0
    while(z<count):
        z+=1
        S+=buildIKF(a,a+h)
        a+=h
    print("s - ",S)
    if k > 3:
        print("m - ", eitkenRule(S, S1, S2, 2))
    E=rungeRule(S1,S,2)
    if E<maxE:
        maxE=E
        maxS = S
    if k>2:

        S2=S1

    S1 = S #первоначальное приближение
print(maxS,maxE,k)
h=2.5
k=0
S = [0,0,0]
count = 1
while(k!=3):
    Sthis = 0
    h = h/2
    a = 0.7
    b = 3.2
    count = count * 2
    z=0
    while(z<count):
        Sthis+=buildIKF(a,a+h)
        a+=h
        z+=1
    S[k]=Sthis
    k += 1 # нахождение hopt
print("done")
m =eitkenRule(S[0], S[1], S[2], 2)
hopt = 0.95*h*2*(10**(-6)/rungeRule(S[1],S[2],2,m))**(1/m)
h = 2.5/round(2.5/hopt)
print(h)
count  = round(2.5/h)
k=0
maxE = 100
while(abs(maxE)>10**-6):

    k+=1
    S = 0
    a = 0.7
    b = 3.2
    z = 1
    while (z < count):
        if (a+h>b):
            print("Smth wrong")
        S+=buildIKF(a,a+h)
        a+=h
        z+=1
    if k > 3:
        m = eitkenRule(S, S1, S2, 2)
        print("m - ", m)
        E=rungeRule(S1,S,2,m)
    if E<maxE:
        maxE=E
        maxS = S
    if k>2:
        S2=S1
    print(E)
    S1 = S
    h = h/2#приближение с hopt
    count = count * 2
print(maxS,maxE,k+3)
a = 0.7
b = 3.2


S1 = buildGaus(a,b)
print(S1)
print(methodicLatencyGaus())
print(analiticLatency(S1))
h=2.5
maxE = 100
k = 1
maxS = 1
count = 1
while(abs(maxE)>10**-6):

    k+=1
    S = 0
    h = h/2
    a = 0.7
    b = 3.2
    count = count * 2
    z = 0
    while (z < count):
        S+=buildGaus(a,a+h)
        a+=h
        z+=1
    if k > 3:
        print("m - ", eitkenRule(S, S1, S2, 2))
    E=rungeRule(S1,S,2)
    if E<maxE:
        maxE=E
        maxS = S
    if k>2:

        S2=S1
    S1 = S
print(maxS,maxE,k)
h=2.5
k=0
S = [0,0,0]
count = 1
while(k!=3):
    Sthis = 0
    h = h/2
    a = 0.7
    b = 3.2
    count = count * 2
    z = 0
    while (z < count):
        Sthis+=buildGaus(a,a+h)
        a+=h
        z+=1
    S[k]=Sthis
    k += 1
m =eitkenRule(S[0], S[1], S[2], 2)
hopt = 0.95*h*2*(10**(-6)/rungeRule(S[1],S[2],2))**(1/m)
h = 2.5/round(2.5/hopt)
print("hopt - ",h)
k=1
maxE = 100
count = round(2.5/hopt)
while(abs(maxE)>10**-6):

    k+=1
    S = 0
    a = 0.7
    b = 3.2
    z = 0
    while (z < count):
        S+=buildGaus(a,a+h)
        a+=h
        z+=1
    if k > 3:
        print("m - ", eitkenRule(S, S1, S2, 2))
    E=rungeRule(S1,S,2)
    if E<maxE:
        maxE=E
        maxS = S
    if k>2:
        S2=S1
    S1 = S
    h = h/2
    count = count * 2
print(maxS,maxE,k+3)