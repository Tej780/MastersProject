import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline as cs
from numba import jit
import warnings
import cProfile
import pyopencl as cl

warnings.simplefilter('ignore', np.RankWarning)


def load(directory):
    data = np.loadtxt(directory, skiprows=1)
    x = directory.split('_')
    y = x[-1].split('x')
    z = np.zeros(len(data[0]))
    z[0] = int(y[0])
    data = np.vstack((data, z))
    return data, z[0]


def NK():
    T = np.arange(0.0, 2.0, 0.01)
    rho = 2 * T / np.pi
    plt.plot(T, rho)


@jit
def xy(data):
    x = []
    y = []
    for i in data[:-1]:
        x.append(i[0])
        y.append(i[1])
    return x, y


def interpolate(x, y, Tbkt, end):
    f = np.poly1d(np.polyfit(x, y, 5))
    xfit = np.linspace(Tbkt, end, 1000)
    yfit = f(xfit)
    return xfit, yfit


@jit
def reduce(upsilon, T, x, y):
    ups = []
    newT = []
    newX = []
    newY = []
    for i in range(len(T)):
        if upsilon[i] > 0:
            ups.append(upsilon[i])
            newT.append(T[i])
            newX.append(x[i + 1])
            newY.append(y[i + 1])
    return ups, newT, newX, newY


def plot(data):
    x, y = xy(data)
    x = x[3:]
    y = y[3:]
    plt.plot(x, y, '--')
    return x, y


@jit
def ups(L, C, c, Tbkt, Tarray):
    Upsilon = []
    Temp = []
    for T in Tarray:
        if (T - Tbkt > 0):

            Ups = (2 / np.pi) * (
                    T / (1 - 0.5 * np.sqrt(c * (T - Tbkt)) / np.tan((np.log(L) + C) * np.sqrt(c * (T - Tbkt)))))
            Temp.append(T)
            Upsilon.append(Ups)

            if (Ups < 0):
                return Upsilon, Temp
    return Upsilon, Temp


@jit
def leastSquare(x, y, ups):
    lsq = 0;
    for i in range(len(x)):
        res = y[i] - ups[i]
        sqr = np.square(res)
        lsq += sqr
    lsq = lsq / len(x)
    return lsq


@jit
def minimize(x, y, L):
    # Arbitrarily pick constants
    print(L)
    Tbkt = 0.8935
    end = 1.5
    Tarray = np.linspace(Tbkt, end, 1000)

    X, Y = interpolate(x, y, Tbkt, end)

    # loop over different constants
    Carray = np.linspace(0.0002, 0.1, 100)
    carray = np.linspace(0.1, 5, 100)

    minlsq = np.inf
    minups = []
    minY = []
    minX = []
    minc = 0
    minC = 0
    for C in Carray:
        for c in carray:
            upsilon = []
            T = []
            X, Y = interpolate(x, y, Tbkt, end)
            # work out Upsilon and T
            upsilon, T = ups(L, C, c, Tbkt, Tarray)

            # Reduce to arrays of same length
            upsilon, T, X, Y = reduce(upsilon, T, X, Y)
            lsq = leastSquare(X, Y, upsilon)
            if (lsq < minlsq):
                minlsq = lsq
                minups = upsilon
                minY = Y
                minX = X
                minC = C
                minc = c
    print(minC)
    print(minc)
    print(minups)
    return minups, minX, minX, minY


def arbitrary(x, y, L):
    # Arbitrarily pick constants
    Tbkt = 0.8935
    end = 1.5
    Tarray = np.linspace(Tbkt, end, 1000)

    # loop over different constants
    C = 0.01
    c = 5
    X, Y = interpolate(x, y, Tbkt, end)
    # work out Upsilon and T
    upsilon, T = ups(L, C, c, Tbkt, Tarray)

    # Reduce to arrays of same length
    upsilon, T, X, Y = reduce(upsilon, T, X, Y)
    lsq = leastSquare(X, Y, upsilon)
    print(lsq)
    return upsilon, T, X, Y


def minFit(directory):
    data, L = load(directory)
    x, y = plot(data)

    # Draw Nelson-Kosterlitz line
    NK()

    upsilon, T, x, y = minimize(x, y, L)
    # upsilon, T, x, y = arbitrary(x, y)

    plt.plot(T, upsilon)
    plt.plot(x, y)

    plt.grid(True)
    plt.show()

def fit(directory):
    data, L = load(directory)
    x, y = plot(data)

    # Draw Nelson-Kosterlitz line
    NK()

    upsilon, T, x, y = arbitrary(x, y, L)
    # upsilon, T, x, y = arbitrary(x, y)

    plt.plot(T, upsilon,'-')
    plt.plot(x, y)

    plt.grid(True)
    plt.show()

'''
##checking differences between 1 run and 2 run datasets
data3 = load("/home/tejan/Desktop/ECMCData/Run5/topology_ECMC_XY_032x032.dat")
data2 = load("/home/tejan/Desktop/ECMCData/Run32-Part2/topology_ECMC_XY_0032x0032.dat")
data1 = load("/home/tejan/Desktop/ECMCData/Run32-Part1/topology_ECMC_XY_0032x0032.dat")

plot(data1)
plot(data2)
plot(data3)

Tbkt = 0.8935
end = 1.5
L = 32
C = 1
c = 3
Tarray = np.linspace(Tbkt, end, 1000)
upsilon, T = ups(L, C, c, Tbkt, Tarray)

plt.plot(T, upsilon)
NK()
plt.grid(True)
plt.show()
'''

##Staring fitting procedure
# load and plot data L=64

#fit("/home/tejan/Desktop/ECMCData/TransitionRegion/256/topology_ECMC_XY_0256x0256.dat")
#minFit("/home/tejan/Desktop/ECMCData/TransitionRegion/256/topology_ECMC_XY_256x256.dat")
fit("/home/tejan/Desktop/ECMCData/TransitionRegion/16/topology_ECMC_XY_016x016.dat")
#minFit("/home/tejan/Desktop/ECMCData/TransitionRegion/16/topology_ECMC_XY_016x016.dat")


