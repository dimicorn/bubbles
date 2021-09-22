import openpyxl
import matplotlib.pyplot as plt
from bubble import*


def data():
    count = 2
    wb = openpyxl.Workbook()
    wb['Sheet'].title = 'Sheet1'
    sh1 = wb.active
    sh1['A1'].value = 'gamma'
    sh1['B1'].value = 'k_rho'
    sh1['C1'].value = 'n_int'
    sh1['D1'].value = 'distance'
    sh1['E1'].value = 'Q_p'

    # This works
    '''
    sh1['D1'].value = 'lambda_c1'
    sh1['E1'].value = 'v1(lambda_c1)'
    sh1['F1'].value = 'v2(lambda_c1)'
    sh1['G1'].value = 'lambda_c2'
    sh1['D1'].value = 'delta lambda %'
    sh1['J1'].value = 'dv, r = R_s'
    sh1['K1'].value = 'dv, r = R_c'
    sh1['L1'].value = 'R_sw (au), v_int = ' + str(v_int) + ', n = ' + str(n0) + ' years, R_c/R_sw = ' + str(k0)
    sh1['F1'].value = 'P(r) from B1'
    sh1['G1'].value = 'P_sw from B2'
    '''

    for i in range(len(G)):
        for j in range(len(K)):
            for q in range(len(N_int)):
                bubble = Bubble(G[i], K[j], N_int[q])
                if bubble.intersect:
                    bubble.values(sh1, count)
                    count += 1
    wb.save(file_path)


def vel_plots(x):
    # Proper code for the velocity plots
    if x == 1:
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
                        x, y = bubble.norm2()
                        plt.plot(x, y, label=Gamma + str(G[i]))
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
            plt.savefig('V(Lambda), ' + N + str(N_int[q]) + '.png', dpi=300)
            plt.show()

    # Searching for intersection code
    elif x == 2:
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    lamb, vel = bubble.lambda_c, bubble.velocity
                    plt.plot(lamb, vel, label=N + str(N_int[q]))
                    if q == 0:
                        plt.plot(bubble.slope_x, bubble.slope_y, label=Gamma + str(G[i]))
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


def pres_plots(x):
    # Repeating plot legend
    def leg():
        plt.plot([], [], 'k-', label=N + str(N_int[0]))
        plt.plot([], [], 'k--', label=N + str(N_int[1]))
        plt.plot([], [], 'k-.', label=N + str(N_int[2]))
        plt.plot([], [], 'k.', label=N + str(N_int[3]))

    # Repeating plot label
    def lab(y):
        plt.grid()
        plt.xlabel('Lambda')
        plt.ylabel('P')
        plt.title('P(Lambda), ' + Gamma + str(G[y]))
        plt.legend()
        plt.savefig('P(Lambda), ' + Gamma + str(G[y]) + '.png', dpi=300)
        plt.show()

    # Code for all the pressure plots
    if x == 1:
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    lamb, pres = bubble.lambda_c, bubble.pressure
                    if q == 0:
                        plt.plot(lamb, pres, color=Colors[j], label=K_rho + str(K[j]), linestyle=Line_style[q])
                    else:
                        plt.plot(lamb, pres, color=Colors[j], linestyle=Line_style[q])
            leg()
            lab(i)

    # Code for the plots that don't have intersection
    elif x == 2:
        for i in range(len(G)):
            for j in range(len(K)):
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    if bubble.intersect:
                        lamb, pres = bubble.lambda_c, bubble.pressure
                        plt.plot(lamb, pres, label=K_rho + str(K[j]) + ', ' + N + str(N_int[q]))
            lab(i)

    # Code for all the plots with intersection
    elif x == 3:
        for i in range(len(G)):
            for j in range(len(K)):
                count = 0
                for q in range(len(N_int)):
                    bubble = Bubble(G[i], K[j], N_int[q])
                    if bubble.intersect:
                        lamb, pres = bubble.lambda_c, bubble.pressure
                        if count == 0:
                            plt.plot(lamb, pres, color=Colors[j], label=K_rho + str(K[j]), linestyle=Line_style[q])
                            count += 1
                        else:
                            plt.plot(lamb, pres, color=Colors[j], linestyle=Line_style[q])
            leg()
            lab(i)


# Creating a scatter plot with values of Q_p in 3D
def scat():
    for i in range(len(G)):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        for j in range(len(K)):
            for q in range(len(N_int)):
                bubble = Bubble(G[i], K[j], N_int[q])
                if bubble.intersect:
                    q_p = bubble.q_p()
                    if 0.9 < q_p < 1.15:
                        ax.scatter3D(K[j], N_int[q], bubble.q_p(), c='blue', alpha=0.5)  # , c=zdata, cmap='Greens'
                    else:
                        ax.scatter3D(K[j], N_int[q], bubble.q_p(), c='red', alpha=0.5)
        plt.title('Q_p, ' + Gamma + str(G[i]))
        ax.set_xlabel('K_rho')
        ax.set_ylabel('N_int')
        ax.set_zlabel('Q_p')
        fig.savefig('Q_p, ' + Gamma + str(G[i]) + '.png', dpi=300)
        fig.show()


# Creating colormaps
def colormap():
    for i in range(len(G)):
        for j in range(len(K)):
            for q in range(len(N_int)):
                bubble = Bubble(G[i], K[j], N_int[q])
                if bubble.intersect:
                    q_p = bubble.q_p()
                    if 0.9 < q_p < 1.15:
                        plt.scatter(K[j], N_int[q], c='blue', alpha=0.5)
                    else:
                        plt.scatter(K[j], N_int[q], c='red', alpha=0.5)
                else:
                    plt.scatter(K[j], N_int[q], c='gray', alpha=0.5)
        plt.title('Colormap, ' + Gamma + str(G[i]))
        plt.xlabel('K_rho')
        plt.ylabel('N_int')
        plt.savefig('Colormap, ' + Gamma + str(G[i]) + '.png', dpi=300)
        plt.show()


if __name__ == '__main__':
    print('You are running the wrong file!')
