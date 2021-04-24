import numpy as np
from Get_global_value import J_type
from Get_global_value import cc
from Get_global_value import Ez


def f_kin_j(RR, AA, q, joints):
    n = len(joints)
    POS_j = np.zeros((n, 3))
    ORI_j = np.zeros((n, 3, 3))
    for i in range(n):
        PorR = (J_type[joints[i]] == 'P')
        ORI_tmp = AA[joints[i], :, :]
        POS_tmp = RR[joints[i], :] + ORI_tmp @ (cc[joints[i], joints[i], :] - (PorR * Ez) * q[joints[i]])

        POS_j[i] = POS_tmp
        ORI_j[i] = ORI_tmp

    return POS_j, ORI_j



