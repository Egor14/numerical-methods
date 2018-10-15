import math

def arctg(x):
    r=1
    f=0
    k=0
    while math.fabs(r)>10**(-6)/3.756:
        r= (-1)**k*x**((-1)*(2*k+1))/(2*k+1)
        f+=r
        k+=1
    return math.pi / 2 * math.copysign(1, x) - f

def sh(x):
    r = 1
    f = 0
    k = 0
    while r > 10 ** (-6) / 3.69/3:
        r = x ** (2 * k + 1) / math.factorial(2 * k + 1)
        f += r
        k += 1
    return f

def sqrt(x):
    w0=1
    w1=2
    while math.fabs(w1-w0)>10**(-6)/4.8:
        c=w1
        w1=0.5*(w0+x/w0)
        w0=c
    return w1

x=0.01
for i in range(10):
    f1=math.sinh(math.sqrt(2*x+0.45))/math.atan(6*x+1)
    f2=sh(sqrt(2*x+0.45))/arctg(6*x+1)
    print('f1= ', '%.10f' % f1, '   f2= ', '%.10f' % f2, '  f1-f2= ', '%.10f' % math.fabs(f1-f2))
    x+=0.005