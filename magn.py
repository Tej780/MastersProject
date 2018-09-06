import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt("/home/tejan/Desktop/Run6/magn_ECMC_XY_064x064.dat", skiprows=1)
data2 = np.loadtxt("/home/tejan/Desktop/Run5/magn_ECMC_XY_032x032.dat", skiprows=1)
data3 = np.loadtxt("/home/tejan/Desktop/Run2/magn_ECMC_XY_004x004.dat", skiprows=1)
global x
global y


def plot(data):
    x = []
    y = []
    for i in data:
        x.append(i[0])
        y.append(i[1])
    plt.plot(x, y)
    x = []
    y = []


plot(data1)
plot(data2)
plot(data3)
T = np.arange(0.0, 2.0, 0.01)
M = 2 * T / np.pi
plt.plot(T, M)
plt.grid(True)
plt.show()
