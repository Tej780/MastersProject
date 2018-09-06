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
    L = data[-1][0]
    x, y = xy(data)
    if (L>10):
        x = x[3:]
        y = y[3:]
    else:
        x = x[1:]
        y = y[1:]
    plt.plot(x, y, '--', label='L=' + str(L),)
    plt.grid()
    plt.axis((0,2.0,0,1.0))
    plt.xlabel('Temperature $T$ (K)')
    plt.ylabel(('$\\rho_{s}$'))
    plt.title('Raw Data for Spin Stiffness $\\rho$ vs Temperature T for increasing system size')
    f = interp1d(x, y, kind=3)
    xfit = np.linspace(x[0], x[-1], 100 * len(x))
    yfit = f(xfit)
    plt.legend()
    plt.savefig('../../../ECMC/data/FirstData.png')
    #plt.plot(xfit, yfit)


# working
dataset = []
Larray = [4,8,16,32,64,128,256,512,1048]
for L in Larray:
    data = load("/home/tejan/ECMC/data/NKF1Data/{:0.0f}x{:0.0f}/topology_ECMC_XY_{:03}x{:03}.dat".format(L,L,int(L),int(L)))
    dataset.append(data)
    plot(data)


cguess = [0.4, 0.8, 1.2, 1.8, 5.2,
          20, 7.0, 9.0]
# working
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
    F = 1 + (g / (2 * np.log(L) + C))
    return F


def sublist(x, y):
    xsub = []
    ysub = []
    for i in range(len(x)):
        if 0.87 < x[i] < 0.91:
            xsub.append(x[i])
            ysub.append(y[i])
    return xsub, ysub


def interpolate(x, y):
    xnew = np.linspace(x[0], x[-1], 1000 * len(x))
    f = interp1d(x, y, kind=7)
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


def C(data1, data2, cguess):
    x1, y1 = xy(data1)
    x2, y2 = xy(data2)
    plt.plot(x1, y1, 'o')
    Cguess = np.linspace(cguess-0.35, cguess+0.35, 70);
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
    plt.plot(xsub, y1new, label='L=16')
    plt.plot(xsub, y2new, label='L=32')
    plt.legend()
    plt.plot(x12, y1new[x12index], 'o')
    T = np.arange(0.0, 2.0, 0.0001)
    rho = 2 * T / np.pi
    plt.plot(T, rho,label='NK relation')
    plt.ylabel("Scaled Spin Stiffness $\\rho_{s}(L)/F_{1}(L)$")
    plt.xlabel("Temperature T (K)")
    plt.title("Intersection method for deducing $C(L_{1},L_{2})$ from the Scaled Spin Stiffness")
    plt.show()
    print(Cmin)
    return Cmin, x12


i = 0
Carray = []
Tarray = []
Larray = []
while i < len(dataset) - 1:
    L = dataset[i][-1][0]
    Cmin, T = C(dataset[i], dataset[i + 1], cguess[i])
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
#
plt.plot(Larray, Carray)
plt.title("A graph showing the relationship between $C(L_{1},L_{2})$ and the inverse system size")
plt.xlabel("System Size $\\frac{1}{L}$")
plt.ylabel("Fitting Parameter $C$")
plt.show()
