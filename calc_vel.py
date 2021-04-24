import numpy as np
from Get_global_value import num_q
from Get_global_value import J_type
from Get_global_value import Ez
from Get_global_value import BB
from Get_global_value import m0
from Get_global_value import m
from Get_global_value import mass
from Get_global_value import inertia0
from Get_global_value import inertia
from Get_global_value import cc
from calc_jr import calc_jr
from calc_jt import calc_jt
from cross import cross
from Get_global_value import c0


def calc_vel(A0, AA, v0, w0, q, qd):
    vv = np.zeros((num_q, 3))
    ww = np.zeros((num_q, 3))
    if num_q == 0:
        print('Single body, there is no link')
    else:
        for i in range(num_q):
            if BB[i] == -1:
                A_I_i = AA[i, :, :]
                if J_type[i] == 'R':
                    ww[i, :] = w0[0:3] + A_I_i @ Ez * qd[i]
                    vv[i, :] = v0[0:3] + cross(w0[0:3], A0 @ c0[i, :]) - cross(ww[i, :], A_I_i @ cc[i, i, :])
                else:
                    ww[i, :] = w0[0:3]
                    vv[i, :] = v0[0:3] \
                               + cross(w0[0:3], np.dot(A0, c0[i, :])) \
                               - cross(ww[i, :], np.dot(A_I_i, cc[i, i, :])) \
                               + cross(ww[i, :], np.dot(np.dot(A_I_i, Ez), q[i])) \
                               + np.dot(np.dot(A_I_i, Ez), qd[i])

            else:
                A_I_BB = AA[BB[i], :, :]
                A_I_i = AA[i, :, :]
                if J_type == 'R':
                    ww[i, :] = ww[BB[i], :] + np.dot(np.dot(A_I_i, Ez), qd[i])
                    vv[i, :] = vv[BB[i], :] \
                               + cross(ww[BB[i], :], np.dot(A_I_BB, cc[BB[i], i, :])) \
                               - cross(ww[i, :], np.dot(A_I_i, cc[i, i, :]))
                else:
                    ww[i, :] = ww[BB[i], :]
                    vv[i, :] = vv[BB[i], :] \
                               + cross(ww[BB[i], :], np.dot(A_I_BB, cc[BB[i], i, :])) \
                               - cross(ww[i, :], np.dot(A_I_i, cc[i, i, :])) \
                               + cross(ww[i, :], np.dot(np.dot(A_I_i, Ez), q[i])) \
                               + np.dot(np.dot(A_I_i, Ez), qd[i])

    return vv, ww























