import numpy
import matplotlib.pyplot as plt
import math

g = 9.80665 #ускорение свободного падения
q = 0.77*50/1000
n = 10
Ux = q*g*(math.sqrt(4*math.pi*0.004*n/g))/0.004
print(Ux) #скорость в центре подъёмной силы
m = 7.91/1000 #масса вертолётика
t = float(input('введите время '))#время
ter = int(input('введите ....')) # 1 - если интересны другие графики
T = 20 #температура воздуха в цельсиях
T = T + 273.15 #перевод в кельвины
P = 750 #давление
p = 0.474*P/T #плотность воздуха
S = 2310/1000000 #плщвдь крыла
Cy = 0.21 #коэффициент подъёмной силы
Cx = 0.2 #коэффициент сопротивления
ay = 0 #ускорение по вертикальной оси
ax = 0 #ускорение по горизонтальной оси
dt = 0.001 #дельта t
y_max = -100 #максимальная высота подъёма
y = 0 #высота во время t
Uy = 0 #скорость по вертикольной оси
Fy = 0 #подъёмная сила
Fx = 0 #сила сопротивления

num = []
num1 = []
num2 = []
num3 = []
num4 = []
num5 = []
num6 = []
tpr = 0
thmax = 0

for i in numpy.arange(0,t,dt): #делим прямую на минимальные отрезки и вычисляем её длину
    Fx = Cx * p * S * ((Ux) ** 2) / 2 # вычисляем силу сопротивления
    Fy = Cy * p * S * ((Ux) ** 2) / 2 # вычисляем подъёмную силу
    ay = Fy/m - g # вычисляем вертикальное ускорение
    ax = -Fx/m # вычисляем горизонтальное ускорение
    Ux = Ux + ax * dt # вычисляем горизонтальную скорость
    Uy = Uy + ay*dt # вычисляем вертикальную скорость
    y = y + Uy*dt # вычислем координату
    num.append(ay)
    num1.append(Uy)
    num2.append(Fy)
    num3.append(Fx)
    num4.append(ax)
    num6.append(Ux)
    num5.append(y)
    tpr+=dt
    if y<=0: #после приземления уже нет смысла считать
        y = 0
        break
    if y_max<y: #находим максимум
        y_max = y
        thmax+=dt

print('конечная координата ',y)
print('максимальная высота полёта ',y_max)
print('время приземления ',tpr)
print('время максисмальной точки ',thmax)
num00 = []
td = len(num)
x = 0
for i in range(td):
    x +=dt
    num00.append(x)
plt.title("Зависимоть координаты от времени ")
plt.grid()
plt.xlabel("время")
plt.ylabel("высота в метрах")
plt.plot(num00, num5)
plt.errorbar(num00, num5, xerr = 0.05, yerr=0.05,fmt = 'b-',ecolor='gray', alpha = 0.08)
plt.show()  # зависимость у от времени

if ter==1: #другие графики
    plt.figure(figsize=(8,8))
    plt.title("Зависимоть ау от времени ")
    plt.grid()
    plt.xlabel("время")
    plt.ylabel("ay м/с^2")
    plt.plot(num00,num)
    plt.errorbar(num00, num, xerr=0.08, yerr=0.2, fmt = 'b-',ecolor='gray', alpha = 0.08)
    plt.show() #зависимость ау от времени

    plt.title("Зависимоть Uy от времени ")
    plt.grid()
    plt.xlabel("время")
    plt.ylabel("Uy в м/с")
    plt.plot(num00,num1)
    plt.errorbar(num00, num1, xerr=0.05, yerr=0.15, fmt='b-', ecolor='gray',alpha = 0.08)
    plt.show() #зависимость Uy от времени

    plt.title("Зависимоть Fy от времени ")
    plt.grid()
    plt.xlabel("время")
    plt.ylabel("Fy в Н")
    plt.plot(num00,num2)
    plt.errorbar(num00, num2, xerr=0.005, yerr=0.008, fmt='b-', ecolor='gray',alpha = 0.08)
    plt.show() #зависимость Fy от времени

    plt.title("Зависимоть Fx от времени ")
    plt.grid()
    plt.xlabel("время")
    plt.ylabel("Fx в Н")
    plt.plot(num00,num3)
    plt.errorbar(num00, num3, xerr=0.005, yerr=0.008, fmt='b-', ecolor='gray',alpha = 0.08)
    plt.show() #зависимость Fx от времени

    plt.title("Зависимоть ax от времени ")
    plt.grid()
    plt.xlabel("время")
    plt.ylabel("ax в м/с^2")
    plt.plot(num00,num4)
    plt.errorbar(num00, num4, xerr=0.08, yerr=0.5, fmt='b-', ecolor='gray',alpha = 0.08)
    plt.show() #зависимость ах от времени

    plt.title("Зависимоть Ux от времени ")
    plt.grid()
    plt.xlabel("время")
    plt.ylabel("Ux в м/с")
    plt.plot(num00,num6)
    plt.errorbar(num00, num6, xerr=0.05, yerr=0.5, fmt='b-', ecolor='gray',alpha = 0.08)
    plt.show() #зависимость Ux от времени
