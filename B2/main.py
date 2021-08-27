import numpy as np
import matplotlib.pyplot as plt
import openpyxl
from constants import*


def R(g, k_rho, n_int, v_int, N, x):  # , l_c
    f_rho = (4 * g / (g + 1) ** 2) ** (1 / (g - 1))
    t_y = 3600 * 24 * 365  # seconds in 1 year
    t = N * t_y  # seconds in N years
    n = (2 + n_int) / (5 - k_rho)
    # R_s = 1
    # R_sw = ((3*(g-1)*n+n_int)/(3*g)*(l_c*R_s)**3/(f_rho*t*v_int)*(g-1)/(g+1))**0.5
    r = (g - 1) / (g + 1) * v_int * f_rho * t * (3 * g / (3 * (g - 1) * n + n_int)) / x ** 3
    return r



s = str(input('Data or plot? '))
if s == 'data':
    file_path = "C:/Users/dmitr/Desktop/MIPT/Лето 21/Dr. Shang/Bubbles/Data/Auto2.xlsx"
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
            for k in range(len(N_int)):
                sh1['A' + str(count)].value = G[i]
                sh1['B' + str(count)].value = K[j]
                sh1['C' + str(count)].value = N_int[k]
                R_sw = R(G[i], K[j], N_int[k], v0, n0, k0)
                sh1['D' + str(count)].value = R_sw / au
                count += 1
    wb.save(file_path)

elif s == 'plot':
    D = {}
    numb = np.linspace(n0, nk, 100)
    V = np.linspace(v0, vk, 100)
    ratio = np.linspace(k0, k_k, 100)
    '''
    for i in range(len(G)):
        D['gamma = ' + str(G[i])] = []
        for j in range(len(numb)):
            R_sw = R(G[i], K[0], N_int[0], v0, numb[j], k0)/au
            D['gamma = ' + str(G[i])].append(R_sw)
        plt.plot(numb, D['gamma = ' + str(G[i])])
    V = np.linspace(v0, vk, 100)
    for i in range(len(G)):
        D['gamma = ' + str(G[i])] = []
        for j in range(len(V)):
            R_sw = R(G[i], K[0], N_int[0], V[j], n0, k0) / au
            D['gamma = ' + str(G[i])].append(R_sw)
        plt.plot(V, D['gamma = ' + str(G[i])])'''
    for i in range(len(G)):
        D['gamma = ' + str(G[i])] = {}
        for j in range(len(K)):
            D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])] = {}
            for k in range(len(N_int)):
                D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])] = {}
                D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])]['numb'] = []
                D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])]['vel'] = []
                D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])]['k'] = []
                for q in range(len(numb)):
                    R_sw = R(G[i], K[j], N_int[k], v0, numb[q], k0) / au
                    D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])]['numb'].append(R_sw)
                for q in range(len(V)):
                    R_sw = R(G[i], K[j], N_int[k], V[q], n0, k0) / au
                    D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])]['vel'].append(R_sw)
                for q in range(len(ratio)):
                    R_sw = R(G[i], K[j], N_int[k], v0, n0, ratio[q]) / au
                    D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])]['k'].append(R_sw)
                # plt.plot(numb, D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])]['numb'])
                # plt.plot(V, D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])]['vel'])
                plt.plot(ratio, D['gamma = ' + str(G[i])]['k_rho = ' + str(K[j])]['n_int = ' + str(N_int[k])]['k'])
    plt.grid()
    plt.show()
