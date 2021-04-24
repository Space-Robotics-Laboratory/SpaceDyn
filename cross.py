import numpy as np


def cross(u, v):
    n = np.zeros(3)
    n[0] = u[1] * v[2] - u[2] * v[1]
    n[1] = u[2] * v[0] - u[0] * v[2]
    n[2] = u[0] * v[1] - u[1] * v[0]

    return n
