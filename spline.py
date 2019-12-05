import numpy as np
from numpy import linalg as LA


def h_vector(xn):
    hn = []
    for i in range(len(xn) - 1):
        h = xn[i + 1] - xn[i]
        hn.append(h)
    return hn


def matrix_a(xn, hn):
    A = np.zeros((len(xn), len(xn)), dtype='float')
    A[0][0] = 2 * hn[0]
    for i in range(len(hn) - 1):
        for j in range(len(hn) - 1):
            if j == i:
                A[i+1][j+1] = 2 * (hn[i] + hn[i + 1])
                A[i][j + 1] = hn[i]
                A[i + 1][j] = hn[i]
    index = len(hn) - 2
    A[A.shape[0] - 1][A.shape[1] - 1] = 2 * hn[len(hn) - 1]
    A[A.shape[0] - 2][A.shape[1] - 1] = hn[len(hn) - 1]
    A[A.shape[0] - 1][A.shape[1] - 2] = hn[len(hn) - 1]
    A[0][0] = 2 * hn[0]
    return A

def matrix_v(hn, fxn, dx):
    v = np.zeros(len(xn), dtype='float')
    for i in range(1,len(v)-1):
        v[i] = 3 * ((1/hn[i]) * (fxn[i + 1] - fxn[i]) - (1/hn[i-1]) * (fxn[i] - fxn[i - 1]))
    v[0] = 3 * ((1/hn[0]) * (fxn[1] - fxn[0]) - dx[0])
    v[len(v) - 1] = 3 * (dx[1] - (1/hn[len(hn) - 1] * (fxn[len(fxn) - 1] - fxn[len(fxn) - 2])))
    return v

def matrix_c(mat_a, mat_v):
    c = LA.solve(mat_a,mat_v)
    return c

def matrix_bd(fxn, hn, c):
    b = np.zeros(len(c) - 1, dtype='float')
    d = np.zeros(len(c) - 1, dtype='float')

    shape = b.shape[0] - 1

    for i in range(len(b)):
        b[i] = ((1/hn[i]) * (fxn[i + 1] - fxn[i])) - ((hn[i]/3) * ((2 * c[i]) + c[i + 1]))
    for i in range(len(d)):
        d[i] = (1/(3 * hn[i])) * (c[i + 1] - c[i])
    return b, d

