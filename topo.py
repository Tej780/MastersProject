import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt("/home/tejan/Desktop/ECMCData/Run6/topology_ECMC_XY_064x064.dat", skiprows=1)
data2 = np.loadtxt("/home/tejan/Desktop/ECMCData/TransitionRegion/64/topology_ECMC_XY_064x064.dat", skiprows=1)
data3 = np.loadtxt("/home/tejan/Desktop/ECMCData/TransitionRegion/512/topology_ECMC_XY_0512x0512.dat", skiprows=1)
#data3 = np.loadtxt("/home/tejan/Desktop/Run2/topology_ECMC_XY_004x004.dat", skiprows=1)
global x
global y


def plot(data):
    x = []
    y = []
    for i in data:
        x.append(i[0])
        y.append(i[1])
    plt.plot(x, y,'o')
    x = []
    y = []


plot(data1)
plot(data2)
plot(data3)
T = np.arange(0.0, 2.0, 0.01)
rho = 2 * T / np.pi
plt.plot(T, rho)

plt.grid(True)
plt.show()
