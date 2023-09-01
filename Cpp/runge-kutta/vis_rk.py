import matplotlib.pyplot as plt
import sys


def FigTitle(filename: str)->str:
    # filename format example: gamma_2_k_rho_2_n_int_1.txt
    temp = filename.split('_')
    g = r"$\gamma=$" + temp[1]
    k = r"$k_{\rho}=$" + temp[4]
    n = r"$\eta_{int}=$" + (temp[-1].split('.'))[0]
    title = ", ".join((g, k, n))
    return title


lambd, velocity, pressure  = [], [], []
density, curve_value = [], []
args = sys.argv[1:]
dir = "data/"
file = open(dir + args[0], "r")

for line in file:
    line = line.split(" ")
    lambd.append(float(line[0]))
    velocity.append(float(line[1]))
    pressure.append(float(line[2]))
    density.append(float(line[3]))
    curve_value.append(float(line[4]))

if args[1] == "vel":
    plt.plot(lambd, velocity)
    plt.plot(lambd, curve_value)
    plt.ylabel(r"$V$")
elif args[1] == "pres":
    plt.plot(lambd, pressure)
    plt.ylabel(r"$P$")
elif args[1] == "den":
    plt.plot(lambd, density)
    plt.ylabel(r"$\rho$")

plt.title(FigTitle(args[0]))
plt.xlabel(r"$\lambda$")
plt.grid(True)
plt.savefig("figures/" + args[0][:-3], dpi=300)