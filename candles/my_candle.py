import matplotlib.pyplot as plt
import numpy
import math
import numpy as np
from scipy.fft import fft, fftfreq

def standardize(data):
    for a in range(2):
        span = max(data[0][a]) - min(data[0][a])
        min_ = min(data[0][a])
        for idx in range(len(data)):
            standardize = (max(data[idx][a]) - min(data[idx][a]))/span
            data[idx][a] = [i/standardize + min_ - min([i/standardize
                            for i in data[idx][a]]) for i in data[idx][a]]
    return data



# противофаза
# ui = 1.297838702869751 # относительная температура для первой свечи
# uj = 0.6100276613214819 # относительная температура для второй всечи
# vi = 0.1735207098559884 # относительная концентрация кислорода, поступающего в пламя, для первой свечи
# vj = 0.13266180927595078 # относительная концентрация кислорода, поступающего в пламя, для второй свечи
#  синфазное колебание
ui = 1.297838702869751 # относительная температура для первой свечи
uj = 1.297838702869751 # относительная температура для второй всечи
vi = 0.1735207098559884 # относительная концентрация кислорода, поступающего в пламя, для первой свечи
vj = 0.1735207098559884 # относительная концентрация кислорода, поступающего в пламя, для второй свечи

t = 2 # время
dt = 0.001 # шаг времени

au = 2.7 #
av = 1 #
e = 10**(-3) #
o = 1 #
c = 5 #
mu0 = 0.5 # коэффициент взаимодействия

num = [] # создание нужных массивов
num1 = []
num2 = []
num3 = []
dui = 0 # создание нужных переменных
duj = 0
dvi = 0
dvj = 0
ui_max = 0 # максимальная температура палмени для первой свечи
t_period_kolebanii_i = 0 # период колебаний для первой свечи
uj_max = 0 # максимальная температура пламени для второй свечи
t_period_kolebanii_j = 0 # период колебаний для второй свечи
ion=0.

int=-1.

for i in numpy.arange(0, t, dt): # вычисление дифференциального уравнения
    dui = (1 / e * (-ui + au * vi * math.exp(ui / (1 + ui / c))) - o * (1 + ui / c) + ion*o * mu0 * (1 + uj / c)**4 + int/mu0 * (1 + uj / c)**4) * dt
    duj = (1 / e * (-uj + au * vj * math.exp(uj / (1 + uj / c))) - o * (1 + uj / c) + ion*o * mu0 * (1 + ui / c)**4 + int/mu0 * (1 + uj / c)**4) * dt
    dvi = (1 - vi - av * vi * math.exp(ui / (1 + ui / c))) * dt
    dvj = (1 - vj - av * vj * math.exp(uj / (1 + uj / c))) * dt
    ui += dui
    uj += duj
    vi += dvi
    vj += dvj
    num.append(ui) # добавление переменных в массивы
    num1.append(uj)
    num2.append(vi)
    num3.append(vj)

dlen = len(num)
time = []
x = 0
for i in range(dlen): # создане массива времени
    x += dt
    time.append(x)

data = [[time,num], [time, num1]]
limits = [(min(data[1][a]), max(data[1][a])) for a in range(2)]
norm_data = standardize(data)

fig, ax = plt.subplots()

for x, y in norm_data:
    ax.plot(x, y)

ax2, ax3 = ax.twinx(), ax.twiny()
ax2.set_ylim(limits[1])
ax3.set_xlim(limits[0])
plt.title('зависимость относительной температуры от времени для двух свечей ')
plt.show()

data = [[time,num2], [time, num3]]
limits = [(min(data[1][a]), max(data[1][a])) for a in range(2)]
norm_data = standardize(data)

fig, ax = plt.subplots()

for x, y in norm_data:
    ax.plot(x, y)

ax2, ax3 = ax.twinx(), ax.twiny()
ax2.set_ylim(limits[1])
ax3.set_xlim(limits[0])

plt.title('зависимость относительной концентрации кислорода для двух всечей ')
plt.show()

print(f'{num[0] = }')
print(f'{num1[100] = }')
print(f'{num2[0] = }')
print(f'{num3[100] = }')

# Number of sample points
N = int(t/dt)
# sample spacing
T = dt

yf = fft(num)
xf = fftfreq(N, T)[:N//2]
plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
plt.grid()
plt.show()