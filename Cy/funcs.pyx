import numpy as np
from scipy import integrate


G = np.array([1.001, 1.01, 1.4, 1.5, 1.6, 5/3, 1.7, 1.9, 2, 4])
K = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3])
N_int = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3])
cdef double x0 = 1
cdef double x_k = 0.7
cdef int v0 = 1
cdef int rho0 = 1
cdef int p0 = 1

class Bubble(object):
    def __init__(self, double a, b, c):
        self.gamma = a
        self.k_rho = b
        self.n_int = c
        cdef double g = self.gamma
        cdef double k = self.k_rho

        # Solving system of ODE's
        def solver(double x, r):
            v, rho, p = r
            cdef double n1 = self.n_int
            cdef double n = (2 + n1) / (5 - k)
            cdef double fv = (4 * g * v / ((g + 1) * x) - k - (1 - 1 / n) * (
                    ((g + 1) / (g - 1)) * (2 * v / (g + 1) - x) * rho * v / p - 2)) * (
                         ((g + 1) / (g - 1)) * (2 * v / (g + 1) - x) ** 2 * rho / p - 2 * g / (g + 1)) ** (-1)
            cdef double fp = ((1 / n - 1) * v - (2 * v / (g + 1) - x) * fv) * ((g + 1) / ((g - 1) * rho))
            cdef double f_rho = ((-k * (g - 1) - 2 * (1 - 1 / n)) * (2 * v / (g + 1) - x) ** (-1) - 1 / p * fp) * (-rho / g)
            return fv, f_rho, fp
        sol = integrate.solve_ivp(solver, (x0, x_k), (v0, rho0, p0),
                                  t_eval=np.linspace(x0, x_k, 1000), method='Radau')
        self.velocity, self.density, self.pressure = np.array(sol.y)
        self.lambda_c = np.array(sol.t)
        self.slope_x = np.linspace(x0, x_k, 1000)
        self.slope_y = np.array([(self.gamma + 1) / 2 * self.slope_x[e] for e in range(1000)])
        self.dist = -1
        self.intersect = False
        cdef int count = 0
        x1, y1, x2, y2 = self.lambda_c, self.velocity, self.slope_x, self.slope_y
        cdef int i = 0
        cdef int j = 0
        cdef int q = 0
        cdef double lamb, vel, s, ox, oy, d
        for i in range(len(x1)):
            for j in range(len(x2)):
                if x1[i] == x2[j] and y1[i] == y2[j]:
                    count += 1
                    self.dist = 0
        if count == 0:
            lamb, vel = self.lambda_c[-1], self.velocity[-1]
            s = (lamb - 2 / (g + 1) * vel) * (2 * (g + 1)) / (4 + (g + 1) ** 2)
            ox = lamb - s * (g + 1) / 2
            oy = ox * (g + 1) / 2
            d = np.sqrt((ox - lamb) ** 2 + (oy - vel) ** 2)
            self.dist = d
        if self.dist < 0.001:
            self.intersect = True
    def q_p(self):
        cdef double g = self.gamma
        cdef double p_r = 2 / (g + 1) * self.pressure[-1]
        cdef double p_sw = (2/(g + 1)) * ((g + 1)**2/(4*g))**(g/(g-1))
        return p_sw/p_r

def loop():
    cdef int i = 0
    cdef int j = 0
    cdef int q = 0
    cdef int len_g = len(G)
    cdef int len_k = len(K)
    cdef int len_n = len(N_int)
    intersects = np.full((len_g, len_k, len_n), -1)
    qp_value = np.full((len_g, len_k, len_n), -1.0)
    for i in range(len_g):
        for j in range(len_k):
            for q in range(len_n):
                bubble = Bubble(G[i], K[j], N_int[q])
                intersects[i][j][q] = bubble.intersect
                qp_value[i][j][q] = bubble.q_p()
    np.save('intersects.npy', intersects)
    np.save('Qp_value.npy', qp_value)