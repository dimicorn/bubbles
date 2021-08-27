from constants import*


class Bubble(object):
    def __init__(self, x, y, z):
        self.gamma = x
        self.k_rho = y
        self.n_int = z

    def r_sw(self, v_int, numb, x):  # , l_c
        g = self.gamma
        n_int = self.n_int
        k_rho = self.k_rho
        f_rho = (4 * g / (g + 1) ** 2) ** (1 / (g - 1))
        t = numb * t_y  # seconds in N years
        n = (2 + n_int) / (5 - k_rho)
        # R_s = 1
        # R_sw = ((3*(g-1)*n+n_int)/(3*g)*(l_c*R_s)**3/(f_rho*t*v_int)*(g-1)/(g+1))**0.5
        r = (g - 1) / (g + 1) * v_int * f_rho * t * (3 * g / (3 * (g - 1) * n + n_int)) / x ** 3
        return r / au


if __name__ == '__main__':
    print('You are running the wrong file!')
