import numpy as np
import matplotlib.pyplot as plt
import openpyxl
from objects import*


def main():
    s = str(input('Data or plot? '))
    if s == 'data':
        wb = openpyxl.Workbook()
        wb['Sheet'].title = "Sheet1"
        sh1 = wb.active
        sh1['A1'].value = 'gamma'
        sh1['B1'].value = 'k_rho'
        sh1['C1'].value = 'n_int'
        sh1['D1'].value = 'R_sw (au)'
        count = 2
        # numb = np.linspace(n0, nk, 100)  # numb of years
        # V = np.linspace(v0, vk, 100)  # velocity
        # K = np.linspace(k0, k_k, 100)
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    sh1['A' + str(count)].value = G[i]
                    sh1['B' + str(count)].value = K[j]
                    sh1['C' + str(count)].value = N_int[q]
                    bubble = Bubble(G[i], K[j], N_int[q])
                    sh1['D' + str(count)].value = bubble.r_sw(v0, numb0, k0)
                    count += 1
        wb.save(file_path)
    elif s == 'plot':
        numb = np.linspace(numb0, numb_k, 100)
        v = np.linspace(v0, vk, 100)
        k = np.linspace(k0, k_k, 100)
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    r_n = []
                    r_v = []
                    r_k = []
                    for w in range(len(numb)):
                        r = bubble.r_sw(v0, numb[w], k0)
                        r_n.append(r)
                    for w in range(len(v)):
                        r = bubble.r_sw(v[w], numb0, k0)
                        r_v.append(r)
                    for w in range(len(k)):
                        r = bubble.r_sw(v0, numb0, k[w])
                        r_k.append(r)
                    plt.plot(numb, r_n)
                    plt.plot(v, r_v)
                    plt.plot(k, r_k)
        plt.grid()
        plt.show()


if __name__ == '__main__':
    main()
