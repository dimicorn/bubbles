import matplotlib.pyplot as plt
import sys
import os


def FigTitle(filename: str) -> str:
    # filename format example: gamma_2_k_rho_2_n_int_1.txt
    temp = filename.split('_')
    g = r"$\gamma=$" + temp[1]
    k = r"$k_{\rho}=$" + temp[4]
    n = r"$\eta_{int}=$" + (temp[-1].split('.'))[0]
    title = ", ".join((g, k, n))
    return title

def main() -> None:
    lambd, velocity, pressure  = [], [], []
    density, curve_value = [], []
    args = sys.argv[1:]

    if len(args) == 0:
        print(f"Add filename with data and choose what you want to plot (\"vel\", \"pres\" or \"den\")!\n"
              f"For example: python3 {sys.argv[0]} gamma_2_k_rho_2_n_int_1.txt vel")
    elif len(args) == 1:
        print(f"Check that you've added filename with data AND have chosen what you want to plot (\"vel\", \"pres\" or \"den\") as command line arguments!\n"
              f"For example: python3 {sys.argv[0]} gamma_2_k_rho_2_n_int_1.txt vel")
    else:
        dir = "data"
        file = open(f"{dir}/{args[0]}", "r")

        output_dir = "figures"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

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
        plt.savefig(f"{output_dir}/{args[0][:-3]}", dpi=300)

if __name__ == "__main__":
    main()
