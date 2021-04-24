import numpy as np
from math import cos
from math import sin
from Get_global_value import d_time


# aw( w0 ) returns a 3x3 transformation representing a rotation about the vector w0.
def aw(omega0):
    omega0_norm = np.linalg.norm(omega0)
    transform_about_omega0 = np.zeros((3, 3))
    if omega0_norm == 0:
        transform_about_omega0 = np.eye(3)
    else:
        theta = omega0_norm * d_time
        unit_omega0 = np.true_divide(omega0, omega0_norm)
        transform_about_omega0[0, 0] = cos(theta) + unit_omega0[0] ** 2 * (1-cos(theta))
        transform_about_omega0[0, 1] = unit_omega0[0] * unit_omega0[1] * (1-cos(theta)) - unit_omega0[2] * sin(theta)
        transform_about_omega0[0, 2] = unit_omega0[2] * unit_omega0[0] * (1-cos(theta)) + unit_omega0[1] * sin(theta)
        transform_about_omega0[1, 0] = unit_omega0[0] * unit_omega0[1] * (1-cos(theta)) + unit_omega0[2] * sin(theta)
        transform_about_omega0[1, 1] = cos(theta) + unit_omega0[1] ** 2 * (1-cos(theta))
        transform_about_omega0[1, 2] = unit_omega0[2] * unit_omega0[1] * (1-cos(theta)) - unit_omega0[0] * sin(theta)
        transform_about_omega0[2, 0] = unit_omega0[2] * unit_omega0[0] * (1-cos(theta)) - unit_omega0[1] * sin(theta)
        transform_about_omega0[2, 1] = unit_omega0[2] * unit_omega0[1] * (1-cos(theta)) + unit_omega0[0] * sin(theta)
        transform_about_omega0[2, 2] = cos(theta) + unit_omega0[2] ** 2 * (1-cos(theta))

    return transform_about_omega0









