import openpyxl
import matplotlib.pyplot as plt
from objects import*
from constants import*


def main():
    a = str(input('Data or plot? '))
    if a == 'plot':
        b = str(input('Velocity or pressure? '))

    if a == 'data':
        count = 2
        wb = openpyxl.Workbook()
        wb['Sheet'].title = 'Sheet1'
        sh1 = wb.active
        sh1['A1'].value = 'gamma'
        sh1['B1'].value = 'k_rho'
        sh1['C1'].value = 'n_int'
        sh1['D1'].value = 'lambda_c1'
        sh1['E1'].value = 'v1(lambda_c1)'
        sh1['F1'].value = 'v2(lambda_c1)'
        sh1['G1'].value = 'lambda_c2'
        sh1['H1'].value = 'delta lambda %'
        sh1['I1'].value = 'dv, r = R_s'
        sh1['J1'].value = 'dv, r = R_c'
        sh1['K1'].value = 'R_sw (au), v_int = ' + str(v_int) + ', n = ' + str(n0) + ' years, R_c/R_sw = ' + str(k0)
        sh1['L1'].value = 'p(lambda_c1)'
        sh1['M1'].value = 'p_sw * rho * v**2'
        sh1['N1'].value = 'f_Psa'
        sh1['O1'].value = '4'+'$/pi /xi$'
        sh1['P1'].value = 'f_p'
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    bubble.values(sh1, count)
                    count += 1
        wb.save(file_path)

    elif a == 'plot' and b == 'vel':
        for q in range(len(N_int)):
            for j in range(len(K)):
                for i in range(len(G)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    lamb, vel = bubble.norm1()
                    if i == 0:
                        plt.plot(lamb, vel, color=Colors[j], label=K_rho + str(K[j]))
                    else:
                        plt.plot(lamb, vel, color=Colors[j])
                    if j == 0:
                        a = [(G[i] + 1) / 2 * lamb[e] for e in range(len(lamb))]
                        lamb, a = norm2(G[i], lamb, a)
                        plt.plot(lamb, a, label=Gamma + str(G[i]))
            plt.grid()
            plt.legend()
            if N_int[q] == 0.5:
                plt.axis([0, 1., 0.5, 1.3])
            elif N_int[q] == 1:
                plt.axis([0.4, 1., 0.9, 1.25])
            elif N_int[q] == 1.5 or N_int[q] == 2.5:
                plt.axis([0.55, 1., 1, 1.5])
            plt.xlabel('Lambda')
            plt.ylabel('V')
            plt.title('V(Lambda), ' + N + str(N_int[q]))
            plt.show()

    elif a == 'plot' and b == 'pres':
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    lamb, pres = bubble.solution2()
                    if q == 0:
                        plt.plot(lamb, pres, color=Colors[j], label=K_rho + str(K[j]))
                    else:
                        plt.plot(lamb, pres, color=Colors[j])
            plt.grid()
            plt.xlabel('Lambda')
            plt.ylabel('P')
            plt.title('P(Lambda), ' + Gamma + str(G[i]))
            plt.legend()
            plt.show()


if __name__ == '__main__':
    main()
