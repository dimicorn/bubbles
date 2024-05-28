import numpy as np
from scipy import integrate
# cimport cython
# from cython.parallel cimport prange


G = np.array([1.001, 1.01, 1.4, 1.5, 1.6, 5./3, 1.7, 1.9, 2, 4])
# K = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3])
# N_int = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3])
K = np.linspace(0, 3, num=20)
N_int = np.linspace(0, 3, num=20)
cdef int x0 = 1
cdef float x_k = 0.7
cdef int v0 = 1
cdef int rho0 = 1
cdef int p0 = 1

# @cython.final
cdef class calc:
    cdef public float gamma, k_rho, n_int, q_p
    cdef public int intersect

    def __cinit__(self, float a, float b, float c):
        self.gamma = a
        self.k_rho = b
        self.n_int = c
        cdef float g = self.gamma
        cdef float k = self.k_rho
        cdef int count = 0
        cdef int i = 0
        cdef int j = 0
        cdef int q = 0
        cdef int e
        cdef float lamb, vel, s, ox, oy

        # Solving system of ODE's
        def solver(x, r):
            cdef v = r[0]
            cdef rho = r[1]
            cdef p = r[2]
            cdef float n1 = self.n_int
            cdef float n = (2 + n1) / (5 - k)
            cdef float fv = (4 * g * v / ((g + 1) * x) - k - (1 - 1 / n) * (
                    ((g + 1) / (g - 1)) * (2 * v / (g + 1) - x) * rho * v / p - 2)) / (
                    ((g + 1) / (g - 1)) * (2 * v / (g + 1) - x)* (2 * v / (g + 1) - x) * rho / p - 2 * g / (g + 1))
            cdef float fp = ((1 / n - 1) * v - (2 * v / (g + 1) - x) * fv) * ((g + 1) / ((g - 1) * rho))
            cdef float f_rho = ((-k * (g - 1) - 2 * (1 - 1 / n)) / (2 * v / (g + 1) - x) - 1 / p * fp) * (-rho / g)
            return fv, f_rho, fp
        cdef sol = integrate.solve_ivp(solver, (x0, x_k), (v0, rho0, p0),
                                  t_eval=np.linspace(x0, x_k, 1000), method='Radau')
        cdef velocity = sol.y[0]
        cdef density = sol.y[1]
        cdef pressure = sol.y[2]
        cdef lambda_c = sol.t
        cdef slope_x = np.linspace(x0, x_k, 1000)
        cdef slope_y = np.array([(g + 1) / 2 * slope_x[e] for e in np.arange(1000)])
        cdef float dist = -1
        self.intersect = 0
        cdef int len_x1 = len(lambda_c)
        cdef int len_x2 = len(slope_x)
        for i in np.arange(len_x1):
            for j in np.arange(len_x2):
                if lambda_c[i] == slope_x[j] and velocity[i] == slope_y[j]:
                    count += 1
                    dist = 0
                    break
        if count == 0:
            lamb = lambda_c[-1]
            vel = velocity[-1]
            s = (lamb - 2 / (g + 1) * vel) * (2 * (g + 1)) / (4 + (g + 1) * (g + 1))
            ox = lamb - s * (g + 1) / 2
            oy = ox * (g + 1) / 2
            dist = (ox - lamb) * (ox - lamb) + (oy - vel) * (oy - vel)
        if dist < 0.00001:
            self.intersect += 1
        cdef float p_r = 2 / (g + 1) * pressure[-1]
        cdef float p_sw = (2/(g + 1)) * ((g + 1)*(g + 1)/(4*g))**(g/(g-1))
        self.q_p = p_sw/p_r


def loop():
    cdef int i = 0
    cdef int j = 0
    cdef int q = 0
    cdef int len_g = len(G)
    cdef int len_k = len(K)
    cdef int len_n = len(N_int)
    intersects = np.full((len_g, len_k, len_n), -1)
    qp_value = np.full((len_g, len_k, len_n), -1.0)
    for i in np.arange(len_g):
        for j in np.arange(len_k):
            for q in np.arange(len_n):
                    bubble = calc(G[i], K[j], N_int[q])
                    intersects[i][j][q] = bubble.intersect
                    qp_value[i][j][q] = bubble.q_p
    np.save('intersects.npy', intersects)
    np.save('Qp_value.npy', qp_value)