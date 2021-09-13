from scipy import integrate
import numpy as np
from constants import*


class Bubble(object):
    def __init__(self, a, b, c):
        self.gamma = a
        self.k_rho = b
        self.n_int = c

        # Solving system of ODE's
        def solver(x, r):
            v, rho, p = r
            g = self.gamma
            k = self.k_rho
            n1 = self.n_int
            n = (2 + n1) / (5 - k)
            fv = (4 * g * v / ((g + 1) * x) - k - (1 - 1 / n) * (
                    ((g + 1) / (g - 1)) * (2 * v / (g + 1) - x) * rho * v / p - 2)) * (
                         ((g + 1) / (g - 1)) * (2 * v / (g + 1) - x) ** 2 * rho / p - 2 * g / (g + 1)) ** (-1)
            fp = ((1 / n - 1) * v - (2 * v / (g + 1) - x) * fv) * ((g + 1) / ((g - 1) * rho))
            f_rho = ((-k * (g - 1) - 2 * (1 - 1 / n)) * (2 * v / (g + 1) - x) ** (-1) - 1 / p * fp) * (-rho / g)
            return fv, f_rho, fp

        # Solutions for the system
        sol = integrate.solve_ivp(solver, (x0, x_k), (v0, rho0, p0),
                                  t_eval=np.linspace(x0, x_k, 1000), method='Radau')
        self.velocity, self.density, self.pressure = sol.y
        self.lambda_c = sol.t

    # Function for eta
    def eta(self):
        k_rho = self.k_rho
        n_int = self.n_int
        x = (2 + n_int) / (5 - k_rho)
        return x

    # Value of the curve at lambda_c
    def v2(self, x, y):
        z = (x + 1) / 2 * y
        return z

    # Approximation using eqn B8a
    def lambda2(self, x, k_rho):
        w = self.eta()
        t = x ** 3 + 12 * x ** 2 + 8 * x + 1 - 0.5 * (x + 1) * (3 * x + 1) * k_rho - (x + 1) * (4 * x + 1) / w
        u = 2 * x ** 3 + 12 * x ** 2 + 7 * x + 1 - 0.5 * (x + 1) * (3 * x + 1) * k_rho - (x + 1) * (4 * x + 1) / w
        return t / u

    # Gradient of velocity, r = R_s
    def dv1(self, x, k_rho):
        w = self.eta()
        y = (-(7 * x + 3) + (x + 1) * k_rho + 3 * (x + 1) / w) / (x + 1)
        return y

    # Gradient of velocity, r = R_c
    def dv2(self, x, k_rho):
        w = self.eta()
        y = (-2 * (x + 1) + k_rho + 2 / w) / (2 * x / (x + 1))
        return y

    # Difference between lambda from the ODE and the approximation using eqn B8a
    def delta(self, x, y):
        z = abs(x - y) / x * 100
        return z

    # Calculating R_sw
    def r_sw(self, v, numb, k):
        g = self.gamma
        f_rho = (4 * g / (g + 1) ** 2) ** (1 / (g - 1))
        t_y = 3600 * 24 * 365
        t = numb * t_y
        n = (2 + self.n_int) / (5 - self.k_rho)
        r = (g - 1) / (g + 1) * v * f_rho * t * (3 * g / (3 * (g - 1) * n + self.n_int)) / k ** 3
        return r / au

    # Adding all the necessary data to the excel table
    def values(self, sh1, count):
        dist = self.d()
        if dist < 0.001:
            sh1['A' + str(count)].value = self.gamma
            sh1['B' + str(count)].value = self.k_rho
            sh1['C' + str(count)].value = self.n_int

            # This works
            '''
            sh1['D' + str(count)].value = self.lambda_c[-1]
            sh1['E' + str(count)].value = self.velocity[-1]
            sh1['F' + str(count)].value = self.v2(self.gamma, self.lambda_c[-1])
            sh1['G' + str(count)].value = l2
            '''
            l2 = self.lambda2(self.gamma, self.k_rho)
            sh1['D' + str(count)].value = self.delta(self.lambda_c[-1], l2)
            sh1['E' + str(count)].value = dist

            # And this too
            '''
            sh1['J' + str(count)].value = self.dv1(self.gamma, self.k_rho)
            sh1['K' + str(count)].value = self.dv2(self.gamma, self.k_rho)
            sh1['L' + str(count)].value = self.r_sw(v_int, n0, k0)
            '''
            sh1['F' + str(count)].value = self.p_r()
            sh1['G' + str(count)].value = self.p_sw()
            sh1['H' + str(count)].value = self.p_sw()/self.p_r()

        # This is broken
        '''
        lamb, pres = self.lambda_c[-1], self.pressure[-1]
        sh1['N' + str(count)].value = pres
        sh1['O' + str(count)].value = self.p_sw()
        sh1['P' + str(count)].value = self.f_psa()
        sh1['Q' + str(count)].value = self.xi()
        sh1['R' + str(count)].value = self.f_p()'''

    # Scaling for the velocity curves
    def norm1(self):
        lamb, vel, pres = self.lambda_c, self.velocity, self.pressure
        for t in range(len(lamb)):
            lamb[t] = 1 - (1 - lamb[t]) / (self.gamma - 1)
        for t in range(len(vel)):
            vel[t] = 1 + (vel[t] - 1) / (self.gamma - 1)
        return lamb, vel

    # Calculating the distance between the last point and the slope
    def d(self):
        g = self.gamma
        x, y = self.lambda_c[-1], self.velocity[-1]
        k = (x - 2/(g+1)*y) * (2*(g+1))/(4+(g+1)**2)
        ox = x - k * (g+1)/2
        oy = ox * (g+1)/2
        r = np.sqrt((ox-x)**2 + (oy-y)**2)
        return r

    # Pressure on the B1 side
    def p_r(self):
        g = self.gamma
        f = 2/(g+1) * self.pressure[-1]
        return f

    # Pressure on the B2 side
    def p_sw(self):
        g = self.gamma
        f = (2/(g + 1)) * ((g + 1)**2/(4*g))**(g/(g-1))
        return f

    '''
    def f_psa(self):
        k = self.k_rho
        g = self.gamma
        n = self.eta()
        f = 2 * ((4 - k)*n - 1)/((g+1)*(3-k)*n)
        return f

    def xi(self):
        g = self.gamma
        n_int = self.n_int
        k = self.k_rho
        n = self.eta()
        f_psa = self.f_psa()
        lambda_c, v = self.lambda_c[-1], self.velocity[-1]
        x = 9*(g - 1)*n_int/((3-k)*n**2*(3*(g-1)*n + n_int)*f_psa*lambda_c**3)
        return x

    def f_p(self):
        k = self.k_rho
        n_int = self.n_int
        f = (4-k)/(3-k) * (n_int/(1+n_int))
        return f
    '''


# Scaling for the linear slope
def norm2(k, a, b):
    for t in range(len(a)):
        a[t] = 1 - (1 - a[t]) / (k - 1)
    for t in range(len(b)):
        b[t] = 1 + (b[t] - 1) / (k - 1)
    return a, b


if __name__ == '__main__':
    print('You are running the wrong file!')
