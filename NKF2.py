##Working
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline as cs


##working
def load(directory):
    data = np.loadtxt(directory, skiprows=1)
    x = directory.split('_')
    y = x[-1].split('x')
    z = np.zeros(len(data[0]))
    z[0] = int(y[0])
    data = np.vstack((data, z))
    return data


##working
def xy(data):
    x = []
    y = []
    for i in data[:-1]:
        x.append(i[0])
        y.append(i[1])
    return x, y


# working
def plot(data):
    x, y = xy(data)
    plt.plot(x, y, 'o')
    f = np.poly1d(np.polyfit(x, y, 10))
    xfit = np.linspace(x[0], x[-1], 100 * len(x))
    yfit = f(xfit)
    plt.plot(xfit, yfit)


# working
data8 = load("/home/tejan/Desktop/ECMCData/512x512/topology_ECMC_XY_512x512.dat")
data7 = load("/home/tejan/Desktop/ECMCData/256x256/topology_ECMC_XY_256x256.dat")
data6 = load("/home/tejan/Desktop/ECMCData/128x128/topology_ECMC_XY_128x128.dat")
data5 = load("/home/tejan/Desktop/ECMCData/Run6/topology_ECMC_XY_064x064.dat")
data4 = load("/home/tejan/Desktop/ECMCData/Run5/topology_ECMC_XY_032x032.dat")
data3 = load("/home/tejan/Desktop/ECMCData/Run4/topology_ECMC_XY_016x016.dat")
data2 = load("/home/tejan/Desktop/ECMCData/Run3/topology_ECMC_XY_008x008.dat")
data1 = load("/home/tejan/Desktop/ECMCData/Run2/topology_ECMC_XY_004x004.dat")
dataset = [data1, data2, data3, data4, data5, data6, data7, data8]
# working
plot(data1)
plot(data2)
plot(data3)
plot(data4)
plot(data5)
plot(data6)
plot(data7)
plot(data8)
T = np.arange(0.0, 2.0, 0.01)
rho = 2 * T / np.pi
plt.plot(T, rho)

plt.grid(True)
plt.show()


# working
def NK(T):
    return (2.0 * T) / np.pi


def F(C, L):
    g = 1.00202783
    F = 1 + (g / (2 * np.log(L) + C + np.log((C / 2.0) + np.log(L))))
    return F


def sublist(x, y):
    xsub = []
    ysub = []
    for i in range(len(x)):
        if 0.8 < x[i] < 9.2:
            xsub.append(x[i])
            ysub.append(y[i])
    return xsub, ysub


def interpolate(x, y):
    xnew = np.linspace(x[0], x[-1], 1000 * len(x))
    f = interp1d(x, y, kind='cubic')
    ynew = f(xnew)
    return xnew, ynew


def intersection(x, y1, y2):
    minx = 0
    xindex = 0
    # working
    for i in range(len(x)):
        if y2[i] - y1[i] < 0:
            minx = x[i]
            xindex = i
            break;
    return minx, xindex


def C(data1, data2):
    x1, y1 = xy(data1)
    x2, y2 = xy(data2)
    plt.plot(x1, y1, 'o')
    Cguess = np.linspace(0.0, 15, 1000);
    Cmin = 0
    mindif = np.inf
    for C in Cguess:
        F1 = F(C, data1[-1][0])
        F2 = F(C, data2[-1][0])
        y1new = y1 / F1
        y2new = y2 / F2
        x2interp, y2interp = interpolate(x2, y2new)
        x1interp, y1interp = interpolate(x1, y1new)
        xsub, y1new = sublist(x1interp, y1interp)
        xsub, y2new = sublist(x2interp, y2interp)
        x12, x12index = intersection(xsub, y1new, y2new)

        if abs(y1new[x12index] - NK(xsub[x12index])) < mindif:
            Cmin = C
            mindif = abs(y1new[x12index] - NK(xsub[x12index]))

    F1min = F(Cmin, data1[-1][0])
    F2min = F(Cmin, data2[-1][0])
    y1min = y1 / F1min
    y2min = y2 / F2min
    x2interp, y2interp = interpolate(x2, y2min)
    x1interp, y1interp = interpolate(x1, y1min)
    xsub, y1new = sublist(x1interp, y1interp)
    xsub, y2new = sublist(x2interp, y2interp)
    x12, x12index = intersection(xsub, y1new, y2new)
    plt.plot(x1interp, y1interp)
    plt.plot(xsub, y1new, label='1sub')
    plt.plot(xsub, y2new, label='2sub')
    plt.legend()
    plt.plot(x12, y1new[x12index], 'o')
    T = np.arange(0.0, 2.0, 0.0001)
    rho = 2 * T / np.pi
    plt.plot(T, rho)
    plt.show()
    print(Cmin)
    return Cmin, x12


i = 0
Carray = []
Tarray = []
Larray = []
while i < len(dataset) - 1:
    Cmin, T = C(dataset[i], dataset[i + 1])
    L = dataset[i][-1][0]
    Carray.append(Cmin)
    Tarray.append(T)
    Larray.append(1.0 / L)
    i += 1;
print(Carray)
print(Tarray)
print(Larray)
plt.plot(Tarray, Carray)
plt.show()
plt.plot(Larray, Tarray)
plt.show()
plt.plot(Larray, Carray)
plt.show()
