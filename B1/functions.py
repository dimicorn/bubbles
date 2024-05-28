import matplotlib.pyplot as plt
from bubble import *

''' data():
    # This works

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

# Proper code for the velocity plots
'''
def vel_plots():
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
'''


# Repeating plot legend for pressure
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


'''
def pres_plots(x):
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
'''

if __name__ == '__main__':
    print('You are running the wrong file!')
