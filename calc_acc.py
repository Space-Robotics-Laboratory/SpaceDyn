import numpy as np
from Get_global_value import num_q
from Get_global_value import J_type
from Get_global_value import Ez
from Get_global_value import BB
from Get_global_value import cc
from cross import cross
from Get_global_value import c0


def calc_acc(A0, AA, w0, ww, vd0, wd0, q, qd, qdd):
    vd = np.zeros((num_q, 3))
    wd = np.zeros((num_q, 3))

    if num_q == 0:
        print('Single body, there is no link')
    else:
        A_I_0 = A0
        for i in range(num_q):
            if BB[i] == -1:
                A_I_i = AA[i, :, :]
                if J_type[i] == 'R':
                    wd[i, :] = wd0[0:3] + cross(ww[i, :], A_I_i @ Ez * qd[i]) + A_I_i @ Ez * qdd[i]
                    vd[i, :] = vd0[0:3] \
                               + cross(wd0[0:3], A_I_0 @ c0[i, :]) \
                               + cross(w0[0:3], cross(w0[0:3], A_I_0 @ c0[i, :])) \
                               - cross(wd[i, :], A_I_i @ cc[i, i, :]) \
                               - cross(ww[i, :], cross(ww[i, :], A_I_i @ cc[i, i, :]))

                else:
                    wd[i, :] = wd0[0:3]
                    vd[i, :] = vd0[0:3] \
                               + cross(wd0[0:3], np.dot(A_I_0, c0[i, :])) \
                               + cross(w0[i, :], cross(w0[0:3], np.dot(A_I_0, c0[i, :]))) \
                               + cross(wd[i, :], np.dot(np.dot(A_I_i, Ez), q[i])) \
                               + cross(ww[i, :], cross(ww[i, :], np.dot(np.dot(A_I_i, Ez), q[i]))) \
                               + 2 * cross(ww[i, :], np.dot(np.dot(A_I_i, Ez), qd[i])) \
                               + np.dot(np.dot(A_I_i, Ez), qdd[i]) \
                               - cross(wd[i, :], np.dot(A_I_i, cc[i, i, :])) \
                               - cross(ww[i, :], cross(ww[i, :], np.dot(A_I_i, cc[i, i, :])))

            else:
                A_I_BB = AA[BB[i], :, :]
                A_I_i = AA[i, :, :]
                if J_type[i] == 'R':
                    wd[i, :] = wd[BB[i], :] + cross(ww[i, :], np.dot(np.dot(A_I_i, Ez), qd[i])) + np.dot(np.dot(A_I_i, Ez), qdd[i])
                    vd[i, :] = vd[BB[i], :] \
                               + cross(wd[BB[i], :], np.dot(A_I_BB, cc[BB[i], i, :])) \
                               + cross(ww[BB[i], :], cross(ww[BB[i], :], np.dot(A_I_BB, cc[BB[i], i, :]))) \
                               - cross(wd[i, :], np.dot(A_I_i, cc[i, i, :])) \
                               - cross(ww[i, :], cross(ww[i, :], np.dot(A_I_i, cc[i, i, :])))
                else:
                    wd[i, :] = wd[BB[i], :]
                    vd[i, :] = vd[BB[i], :] \
                               + cross(wd[BB[i], :], np.dot(A_I_BB, cc[BB[i], i, :])) \
                               + cross(ww[BB[i], :], cross(ww[BB[i], :], np.dot(A_I_BB, cc[BB[i], i, :]))) \
                               + cross(wd[i, :], np.dot(np.dot(A_I_i, Ez), q[i])) \
                               + cross(ww[i, :], cross(ww[i, :], np.dot(np.dot(A_I_i, Ez), q[i]))) \
                               + 2 * cross(ww[i, :], np.dot(np.dot(A_I_i, Ez), qd[i])) \
                               + np.dot(np.dot(A_I_i, Ez), qdd[i]) \
                               - cross(wd[i, :], np.dot(A_I_i, cc[i, i, :])) \
                               - cross(ww[i, :], cross(ww[i, :], np.dot(A_I_i, cc[i, i, :])))

    return vd, wd




















