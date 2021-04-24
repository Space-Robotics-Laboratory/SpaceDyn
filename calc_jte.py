import numpy as np
from Get_global_value import J_type
from Get_global_value import Ez
from Get_global_value import num_q
from f_kin_j import f_kin_j
from f_kin_e import f_kin_e
from cross import cross


def calc_jte(RR, AA, q, joints):
    n = len(joints)
    JJ_te = np.zeros((n, 3))
    if num_q == 0:
        print('Single body, there is no link')
    else:
        POS_j, ORI_j = f_kin_j(RR, AA, q, joints)
        POS_e, ORI_e = f_kin_e(RR, AA, joints)
        for i in range(n):
            A_I_i = AA[joints[i], :, :]
            if J_type[joints[i]] == 'R':
                temp = cross((A_I_i @ Ez), (POS_e - POS_j[i, :]))
                JJ_te[i, :] = temp
            else:
                JJ_te[i, :] = A_I_i @ Ez

    return JJ_te















