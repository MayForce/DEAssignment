import numpy as np
import matplotlib as mplt
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider


#cleans all subplots, counts new values and draws new graphics
def update():
    global slider
    global exactSol
    global euler
    global eulerError
    global impEuler
    global impEulError
    global rungK
    global rungError

    numbOfPoints = int(slider.val)
    step = (xFin - x0) / numbOfPoints

    exactSol.clear()
    exactSol.grid(True)
    euler.clear()
    euler.grid(True)
    eulerError.clear()
    eulerError.grid(True)
    impEuler.clear()
    impEuler.grid(True)
    impEulError.clear()
    impEulError.grid(True)
    rungK.clear()
    rungK.grid(True)
    rungError.clear()
    rungError.grid(True)

    x = np.linspace(x0, xFin, numbOfPoints)
    y = np.zeros([numbOfPoints])
    exactY = np.zeros([numbOfPoints])
    c = np.exp(x0) * (y0 + x0 - 1)
    for i in range(0, numbOfPoints):
        y[i] = 1 - x[i] + c * np.exp(- x[i])
        exactY[i] = 1 - x[i] + c * np.exp(- x[i])
    exactSol.set_ylabel("y")
    exactSol.set_xlabel("x")
    exactSol.plot(x, y, "black")
    exactSol.set_title("exact solution")


    #Euler method
    x = np.linspace(x0, xFin, numbOfPoints)
    y = np.zeros([numbOfPoints])
    y[0] = y0
    for i in range(1, numbOfPoints):
        y[i] = y[i-1] + step * (-y[i-1] - x[i-1])
    euler.set_ylabel("y")
    euler.set_xlabel("x")
    euler.set_title("Euler method")
    euler.plot(x, y, "blue")

    #Euler method errors
    for i in range(0, numbOfPoints):
        y[i] = y[i] - exactY[i]
    eulerError.set_ylabel("y")
    eulerError.set_xlabel("x")
    eulerError.set_title("Euler method errors")
    eulerError.plot(x, y, "blue")


    #improve Euler method
    x = np.linspace(x0, xFin, numbOfPoints)
    y = np.zeros([numbOfPoints])
    y[0] = y0
    for i in range(1, numbOfPoints):
        y[i] = y[i-1] + step * (- (x[i-1] + step/2) - (y[i-1] + step * (- x[i - 1] - y[i - 1]) / 2))
    impEuler.set_ylabel("y")
    impEuler.set_xlabel("x")
    impEuler.set_title("improve Euler method")
    impEuler.plot(x, y, "yellow")

    #Improve Euler method errors
    for i in range(0, numbOfPoints):
        y[i] = y[i] - exactY[i]
    impEulError.set_ylabel("y")
    impEulError.set_xlabel("x")
    impEulError.set_title("improve Euler method errors")
    impEulError.plot(x, y, "yellow")


    #Runge–Kutta method
    x = np.linspace(x0, xFin, numbOfPoints)
    y = np.zeros([numbOfPoints])
    y[0] = y0
    for i in range(1, numbOfPoints):
        k1 = - x[i - 1] - y[i - 1]
        k2 = - (x[i - 1] + step/2) - (y[i - 1] + step * k1 / 2)
        k3 = - (x[i - 1] + step / 2) - (y[i - 1] + step * k2 / 2)
        k4 = - (x[i - 1] + step) - (y[i - 1] + step * k3)
        y[i] = y[i - 1] + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    rungK.set_ylabel("y")
    rungK.set_xlabel("x")
    rungK.set_title("Runge–Kutta methods")
    rungK.plot(x, y, "red")
    #Runge–Kutta method errors
    for i in range(0, numbOfPoints):
        y[i] = y[i] - exactY[i]
    rungError = fig.add_subplot(426)
    rungError.set_ylabel("y")
    rungError.set_xlabel("x")
    rungError.set_title("Runge–Kutta methods errors")
    rungError.plot(x, y, "red")
    plt.draw()

#trecs changes in slider and then updates graphics
def onChangeValue(value):
    update()

print("Enter x0:")
x0 = float(input())
print("Enter y0:")
y0 = float(input())
print("Enter x:")
xFin = float(input())
print("Enter number of points:")
numbOfPoints = int(input())
step = (xFin - x0) / numbOfPoints
fig = plt.figure(figsize=(8, 7))
fig.subplots_adjust(left=0.125, right=0.9, bottom=0.1, top=0.9, wspace=0.3, hspace=1.0)
exactSol = fig.add_subplot(421)
euler = fig.add_subplot(423)
eulerError = fig.add_subplot(422)
impEuler = fig.add_subplot(425)
impEulError = fig.add_subplot(424)
rungK = fig.add_subplot(427)
rungError = fig.add_subplot(426)
ErrorTotal = fig.add_subplot(428)

xT = np.linspace(2, 500, 499)
yTotalE = np.zeros([499])
yTotalImpr = np.zeros([499])
yTotalR = np.zeros([499])
for k in range(2, 500):
    step = (xFin - x0)/k
    x = np.linspace(x0, xFin, k)
    exactY = np.zeros([k])
    c = np.exp(x0) * (y0 + x0 - 1)
    for i in range(0, k):
        exactY[i] = 1 - x[i] + c * np.exp(- x[i])

    x = np.linspace(x0, xFin, k)
    yE = np.zeros([k])
    yE[0] = y0
    yTotalE[0] = abs(yE[0] - exactY[0])
    for i in range(1, k):
        yE[i] = yE[i - 1] + step * (-yE[i - 1] - x[i-1])
        yTotalE[k - 1] = max(yTotalE[k - 1], abs(yE[i] - exactY[i]))

    yImpr = np.zeros([k])
    yImpr[0] = y0
    yTotalE[0] = abs(yImpr[0] - exactY[0])
    for i in range(1, k):
        yImpr[i] = yImpr[i-1] + step * (- (x[i-1] + step/2) - (yImpr[i-1] + step * (- x[i - 1] - yImpr[i - 1]) / 2))
        yTotalImpr[k - 1] = max(yTotalImpr[k - 1], abs(yImpr[i]-exactY[i]))

    yR = np.zeros([k])
    yR[0] = y0
    yTotalE[0] = abs(yR[0] - exactY[0])
    for i in range(1, k):
        k1 = - x[i - 1] - yR[i - 1]
        k2 = - (x[i - 1] + step / 2) - (yR[i - 1] + step * k1 / 2)
        k3 = - (x[i - 1] + step / 2) - (yR[i - 1] + step * k2 / 2)
        k4 = - (x[i - 1] + step) - (yR[i - 1] + step * k3)
        yR[i] = yR[i - 1] + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        yTotalR[k - 1] = max(yTotalR[k - 1], abs(yR[i] - exactY[i]))

ErrorTotal.grid(True)
ErrorTotal.plot(xT, yTotalE, "blue")
ErrorTotal.plot(xT, yTotalImpr, "yellow")
ErrorTotal.plot(xT, yTotalR, "red")
ErrorTotal.set_title("Total error approximation")

axes_slider = plt.axes([0.65, 0.02, 0.22, 0.02])
slider = Slider(axes_slider, label='number of points', valmin=2, valmax=500, valinit=100)
slider.on_changed(onChangeValue)
update()
plt.show()