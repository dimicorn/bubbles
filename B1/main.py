import openpyxl
from functions import *


def main():
    a = str(input('Data or plot? '))
    c = ''
    d = ''
    if a == 'plot':
        b = str(input('Scatter or no? '))
        if b == 'scat':
            d = str(input('3d or 2d? '))
        elif b == 'no':
            c = str(input('Velocity or pressure? '))

    if a == 'data':
        count = 2
        wb = openpyxl.Workbook()
        wb['Sheet'].title = 'Sheet1'
        sh1 = wb.active
        sh1['A1'].value = 'gamma'
        sh1['B1'].value = 'k_rho'
        sh1['C1'].value = 'n_int'
        sh1['D1'].value = 'distance'
        sh1['E1'].value = 'Q_p'

    for i in range(len(G)):
        if d == '3d':
            fig = plt.figure()
            ax = plt.axes(projection='3d')
        for j in range(len(K)):
            if c == 'pres':
                count = 0
            for q in range(len(N_int)):
                bubble = Bubble(G[i], K[j], N_int[q])
                if a == 'data' and bubble.intersect:
                    bubble.values(sh1, count)
                    count += 1
                elif d == '3d' and bubble.intersect:
                    q_p = bubble.q_p()
                    if 0.9 < q_p < 1.15:
                        ax.scatter3D(K[j], N_int[q], bubble.q_p(), c='blue', alpha=0.5)  # , c=zdata, cmap='Greens'
                    else:
                        ax.scatter3D(K[j], N_int[q], bubble.q_p(), c='red', alpha=0.5)
                elif d == '2d':
                    if bubble.intersect:
                        q_p = bubble.q_p()
                        if 0.9 < q_p < 1.15:
                            plt.scatter(K[j], N_int[q], c='blue', alpha=0.5)
                        else:
                            plt.scatter(K[j], N_int[q], c='red', alpha=0.5)
                    else:
                        plt.scatter(K[j], N_int[q], c='gray', alpha=0.5)
                elif c == 'vel':
                    lamb, vel = bubble.lambda_c, bubble.velocity
                    plt.plot(lamb, vel, label=N + str(N_int[q]))
                    if q == 0:
                        plt.plot(bubble.slope_x, bubble.slope_y, label=Gamma + str(G[i]))
                elif c == 'pres' and bubble.intersect:
                    lamb, pres = bubble.lambda_c, bubble.pressure
                    if count == 0:
                        plt.plot(lamb, pres, color=Colors[j], label=K_rho + str(K[j]), linestyle=Line_style[q])
                        count += 1
                    else:
                        plt.plot(lamb, pres, color=Colors[j], linestyle=Line_style[q])
            if c == 'vel':
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

        if d == '3d':
            plt.title('Q_p, ' + Gamma + str(G[i]))
            ax.set_xlabel('K_rho')
            ax.set_ylabel('N_int')
            ax.set_zlabel('Q_p')
            fig.savefig('Q_p, ' + Gamma + str(G[i]) + '.png', dpi=300)
            fig.show()
        elif d == '2d':
            plt.title('Colormap, ' + Gamma + str(G[i]))
            plt.xlabel('K_rho')
            plt.ylabel('N_int')
            plt.savefig('Colormap, ' + Gamma + str(G[i]) + '.png', dpi=300)
            plt.show()
        elif c == 'pres':
            leg()
            lab(i)
    if a == 'data':
        wb.save(file_path)


if __name__ == '__main__':
    main()
