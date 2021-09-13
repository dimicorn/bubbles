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

        # This works
        '''
        sh1['D1'].value = 'lambda_c1'
        sh1['E1'].value = 'v1(lambda_c1)'
        sh1['F1'].value = 'v2(lambda_c1)'
        sh1['G1'].value = 'lambda_c2'
        '''

        sh1['D1'].value = 'delta lambda %'
        sh1['E1'].value = 'distance'

        # And this too
        '''
        sh1['J1'].value = 'dv, r = R_s'
        sh1['K1'].value = 'dv, r = R_c'
        sh1['L1'].value = 'R_sw (au), v_int = ' + str(v_int) + ', n = ' + str(n0) + ' years, R_c/R_sw = ' + str(k0)
        '''

        sh1['F1'].value = 'P(r) from B1'
        sh1['G1'].value = 'P_sw from B2'
        sh1['H1'].value = 'Q_p'

        # These are broken
        '''
        sh1['N1'].value = 'p(lambda_c1)'
        sh1['O1'].value = 'p_sw * rho * v**2'
        sh1['P1'].value = 'f_Psa'
        sh1['Q1'].value = '4'+'$/pi /xi$'
        sh1['R1'].value = 'f_p'''

        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    if bubble.d() < 0.001:
                        bubble.values(sh1, count)
                        count += 1
        wb.save(file_path)

    elif a == 'plot' and b == 'vel':
        # Proper code for the velocity plots
        '''
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
            plt.show()'''

        # Searching for intersection code
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    lamb, vel = bubble.lambda_c, bubble.velocity
                    plt.plot(lamb, vel, label=N + str(N_int[q]))
                    X = np.linspace(x0, 0.7, 1000)
                    a = [(G[i] + 1) / 2 * X[e] for e in range(1000)]
                    if q == 0:
                        plt.plot(X, a, label=Gamma + str(G[i]))
                if i == 0:
                    plt.axis([0.99, 1., 0.99, 1.02])
                elif i == 1:
                    plt.axis([0.95, 1., 0.9, 1.2])
                elif i == 2:
                    plt.axis([0.9, 1., 0.9, 1.2])
                elif i == 3:
                    plt.axis([0.825, 1., 0.9, 1.2])
                elif i == 4:
                    plt.axis([0.75, 1., 0.9, 1.3])
                plt.grid()
                plt.legend()
                plt.xlabel('Lambda')
                plt.ylabel('V')
                plt.title('V(Lambda), ' + Gamma + str(G[i]) + ', ' + K_rho + str(K[j]))
                plt.savefig('V(Lambda), ' + Gamma + str(G[i]) + ', ' + K_rho + str(K[j]) + '.png', dpi=300)
                plt.show()

    elif a == 'plot' and b == 'pres':
        # Code for all the pressure plots
        '''
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    lamb, pres = bubble.lambda_c, bubble.pressure
                    if q == 0:
                        plt.plot(lamb, pres, color=Colors[j], label=K_rho + str(K[j]), linestyle=Line_style[q])
                    else:
                        plt.plot(lamb, pres, color=Colors[j], linestyle=Line_style[q])
            plt.grid()
            plt.plot([], [], 'k-', label=N + str(N_int[0]))
            plt.plot([], [], 'k--', label=N + str(N_int[1]))
            plt.plot([], [], 'k-.', label=N + str(N_int[2]))
            plt.plot([], [], 'k.', label=N + str(N_int[3]))
            plt.xlabel('Lambda')
            plt.ylabel('P')
            plt.title('P(Lambda), ' + Gamma + str(G[i]))
            plt.legend()
            plt.show()'''

        # Code for the plots that don't have intersection
        '''
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    if bubble.d() > 0.001 and not (G[i] == 5/3 and K[j] == 1 and N_int[q] == 1):
                        lamb, pres = bubble.lambda_c, bubble.pressure
                        plt.plot(lamb, pres, label=K_rho + str(K[j]) + ', ' + N + str(N_int[q]))
            plt.grid()
            plt.xlabel('Lambda')
            plt.ylabel('P')
            plt.title('P(Lambda), ' + Gamma + str(G[i]))
            plt.legend()
            plt.show()'''

        # Code for all the plots with intersection
        for i in range(len(G)):
            for j in range(len(K)):
                count = 0
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    if bubble.d() < 0.001 or (G[i] == 5 / 3 and K[j] == 1 and N_int[q] == 1):
                        lamb, pres = bubble.lambda_c, bubble.pressure
                        if count == 0:
                            plt.plot(lamb, pres, color=Colors[j], label=K_rho + str(K[j]), linestyle=Line_style[q])
                            count += 1
                        else:
                            plt.plot(lamb, pres, color=Colors[j], linestyle=Line_style[q])
            plt.grid()
            plt.plot([], [], 'k-', label=N + str(N_int[0]))
            plt.plot([], [], 'k--', label=N + str(N_int[1]))
            plt.plot([], [], 'k-.', label=N + str(N_int[2]))
            plt.plot([], [], 'k.', label=N + str(N_int[3]))
            plt.xlabel('Lambda')
            plt.ylabel('P')
            plt.title('P(Lambda), ' + Gamma + str(G[i]))
            plt.legend()
            plt.savefig('P(Lambda), ' + Gamma + str(G[i]) + '.png', dpi=300)
            plt.show()


if __name__ == '__main__':
    main()
