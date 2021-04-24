from Get_global_value import J_type
from Get_global_value import num_q
from Get_global_value import Ez
import numpy as np


def calc_jre(AA, joints):
    n = len(joints)
    JJ_re = np.zeros((n, 3))
    if num_q == 0:
        print('Single body, there is no link')
    else:
        for i in range(n):
            A_I_i = AA[joints[i], :, :]
            if J_type[joints[i]] == 'R':
                JJ_re[i, :] = A_I_i @ Ez
            else:
                JJ_re[i, :] = np.zeros(3)

    return JJ_re
