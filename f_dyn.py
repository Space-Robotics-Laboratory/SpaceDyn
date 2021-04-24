import numpy as np
from Get_global_value import SE
from Get_global_value import ce
from Get_global_value import num_q
from calc_aa import calc_aa
from calc_pos import calc_pos
from calc_hh import calc_hh
from r_ne import r_ne
from j_num import j_num
from calc_je import calc_je
from skew_sym import skew_sym


def f_dyn(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau):

    AA = calc_aa(A0, q)
    RR = calc_pos(R0, A0, AA, q)
    HH = calc_hh(R0, RR, A0, AA)
    qdd_tempt = np.zeros(num_q)
    vd0_tempt = np.zeros(3)
    wd0_tempt = np.zeros(3)
    Fe_tempt = np.zeros((num_q, 3))
    Te_tempt = np.zeros((num_q, 3))

    Force0 = r_ne(R0, RR, A0, AA, v0, w0, vd0_tempt, wd0_tempt, q, qd, qdd_tempt, Fe_tempt, Te_tempt)

    Force = np.zeros(6 + num_q)
    Force_ex = np.zeros(6 + num_q)
    Fx = np.zeros(3)
    Tx = np.zeros(3)
    taux = np.zeros(num_q)
    Fe_Te = np.zeros((num_q, 6))
    F_ex = np.zeros((num_q, num_q+6))

    Force[0: 3] = F0
    Force[3: 6] = T0

    if num_q == 0:
        print("Single body, there is no link")
    else:
        Force[6: num_q + 6] = tau
        E_3 = np.eye(3)
        O_3 = np.zeros((3, 3))
        num_e = 0
        Me_i = np.zeros((num_q+6, 6))

        for i in range(num_q):
            if SE[i] == 1:
                joints = j_num(num_e)
                tmp = calc_je(RR, AA, q, joints)
                JJ_tx_i = tmp[:, 0:3]
                JJ_rx_i = tmp[:, 3:6]
                num_e = num_e + 1
                A_I_i = AA[i, 0:3, 0:3]
                Re0i = RR[i, :] - R0 + A_I_i @ ce[i, :]
                Me_i[0:3, 0:3] = E_3
                Me_i[0:3, 3:6] = O_3
                Me_i[3:6, 0:3] = skew_sym(Re0i)

                Me_i[3:6, 3:6] = E_3
                Me_i[6:num_q+6, 0:3] = JJ_tx_i
                Me_i[6:num_q+6, 3:6] = JJ_rx_i
                Fe_Te[i, 0:3] = Fe[i, 0:3]
                Fe_Te[i, 3:6] = Te[i, 0:3]
                F_ex[i, :] = Me_i @ Fe_Te[i, 0:6]
        for i in range(num_q):
            Fx += F_ex[i, 0:3]
            Tx += F_ex[i, 3:6]
            taux += F_ex[i, 6:6 + num_q]

    Force_ex[0:3] = Fx
    Force_ex[3:6] = Tx
    Force_ex[6:num_q+6] = taux
    a_Force = Force - Force0 + Force_ex

    Acc = np.linalg.inv(HH) @ a_Force
    vd0 = Acc[0:3]
    wd0 = Acc[3:6]
    qdd = Acc[6:6 + num_q]

    return vd0, wd0, qdd
