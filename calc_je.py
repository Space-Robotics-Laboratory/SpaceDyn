import numpy as np
from calc_jte import calc_jte
from calc_jre import calc_jre
from Get_global_value import num_q


def calc_je(RR, AA, q, joints):
    n = len(joints)
    JJ = np.zeros((n, 6))
    JJ_te = calc_jte(RR, AA, q, joints)
    JJ_re = calc_jre(AA, joints)
    for i in range(n):
        JJ[i, 0:3] = JJ_te[i, 0:3]
        JJ[i, 3:6] = JJ_re[i, 0:3]

    Jacobian = np.zeros((num_q, 6))
    for i in range(n):
        Jacobian[joints[i], :] = JJ[i, 0:6]

    return Jacobian


