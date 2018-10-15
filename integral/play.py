from functions import *

wolf_result = 10.657229068
R1 = 4.360534
R2 = 0.059406

a = 2.5
b = 4.3
alpfa = 2 / 7
beta = 0

h = 4.3 - 2.5
S = 0
count = 1
error = 1
rh1 = 0
mh1 = 0

print('Ответ который должен быть ', wolf_result, '\n')
print("Посчитать интеграл первым или вторым способом?")

# Первый пункт задания

option = int(input())
if option == 1:
    S1 = IKF(a, b)
    print(S1)
else:
    S1 = Gauss(a, b)
    print(S1)

# Второй пункт задания


S2 = 0
iteration = 0
m = 4
while abs(error) > 10e-6:
    a = 2.5
    S = 0
    h = h / 2
    count *= 2
    for i in range(count):
        if option == 1:
            S += IKF(a, a + h)
        else:
            S += Gauss(a, a + h)
        a += h
    print(S)
    if iteration > 0:
        m = eitken(S2, S1, S, 2)
        print('Скорость сходимости ', m)
    error = runge(S, S1, 2, m)
    if iteration == 1:
        rh1 = error
        mh1 = eitken(S2, S1, S, 2)
    S2 = S1
    S1 = S
    iteration += 1

print('\nОконцчательный ответ ', S)

# Третий пункт задания

S = 0
a = 2.5
h = 4.3 - 2.5
h1 = h / 2
h2 = h / 4
hopt = 0.95 * h1 * (10 ** (-6) / abs(rh1)) ** (1 / mh1)
print('\nHopt = ', hopt)
count = round(h / hopt)
hopt = h / count
for i in range(int(count)):
    if option == 1:
        S += IKF(a, a + hopt)
    else:
        S += Gauss(a, a + hopt)
    a += hopt
print(S)
