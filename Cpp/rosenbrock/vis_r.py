import matplotlib.pyplot as plt
import sys


lambd = []
velocity = []
pressure = []
density = []
curve_value= []
args = sys.argv[1:]
dir = "data/"
file = open(dir + args[0], "r")
for i in file:
    j = i.split(" ")
    lambd.append(float(j[0]))
    velocity.append(float(j[1]))
    pressure.append(float(j[2]))
    density.append(float(j[3]))
    curve_value.append(float(j[4]))

plt.plot(lambd, velocity)
# plt.plot(lambd, pressure)
# plt.plot(lambd, density)
plt.plot(lambd, curve_value)
tit = args[0].split('_')
g = tit[0] + "=" + tit[1]
k = tit[2] + "_" + tit[3] + "=" + tit[4]
n = tit[5] + "_" + tit[6] + "=" + (tit[-1].split('.'))[0]
title = ", ".join((g, k, n))
plt.title(title)
plt.savefig("figures/" + args[0][:-3])
    