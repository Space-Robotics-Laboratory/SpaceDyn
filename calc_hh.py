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
from skew_sym import skew_sym


def calc_hh(R0, RR, A0, AA):
    JJ_t = calc_jt(RR, AA)
    JJ_r = calc_jr(AA)

    wE = mass * np.eye(3)

    JJ_tg = np.zeros((num_q, 3))
    HH_w = np.zeros((3, 3))
    HH_wq = np.zeros((num_q, 3))
    HH_q = np.zeros((num_q, num_q))
    HH = np.zeros((num_q+6, num_q+6))
    wr0g = np.zeros(3)

    if num_q == 0:
        HH_w = A0 @ inertia0 @ A0.T
        HH[0:3, 0:3] = wE
        HH[0:3, 3:6] = np.zeros((3, 3))
        HH[3:6, 0:3] = np.zeros((3, 3))
        HH[3:6, 3:6] = HH_w
    else:
        Rm = np.zeros(3)
        for i in range(num_q):
            Rm += m[i] * RR[i, :]

        Rm += m0 * R0
        Rg = Rm / mass
        wr0g = (Rg - R0) * mass

    for i in range(num_q):
        r0i = RR[i, :] - R0
        A_I_i = AA[i, :, :]
        JJ_tg += m[i] * JJ_t[i, :, :]
        inertia_I_i = A_I_i @ inertia[i, :, :] @ A_I_i.T
        HH_w += inertia_I_i + m[i] * (skew_sym(r0i).T @ skew_sym(r0i))
        HH_wq += (inertia_I_i @ JJ_r[i, :, :].T + m[i] * (skew_sym(r0i) @ JJ_t[i, :, :].T)).T
        HH_q += JJ_r[i, :, :] @ inertia_I_i @ JJ_r[i, :, :].T + m[i] * (JJ_t[i, :, :] @ JJ_t[i, :, :].T)

    HH_w += A0 @ inertia0 @ A0.T
    HH[0:3, 0:3] = wE
    HH[0:3, 3:6] = skew_sym(wr0g).T
    HH[0:3, 6:6+num_q] = JJ_tg.T
    HH[3:6, 0:3] = skew_sym(wr0g)
    HH[3:6, 3:6] = HH_w
    HH[3:6, 6:6+num_q] = HH_wq.T
    HH[6:6+num_q, 0:3] = JJ_tg
    HH[6:6+num_q, 3:6] = HH_wq
    HH[6:6+num_q, 6:6+num_q] = HH_q

    return HH


































