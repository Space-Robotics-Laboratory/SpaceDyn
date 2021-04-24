import numpy as np
from Get_global_value import num_q
from Get_global_value import J_type
from Get_global_value import Ez
from Get_global_value import cc
from Get_global_value import BB
from cross import cross


def calc_jt(RR, AA):
    JJ_t = np.zeros((num_q, num_q, 3))
    if num_q == 0:
        print('Single body, there is no link')
    else:
        for i in range(num_q):
            j = i
            while j > -1:
                A_I_j = AA[j, :, :]
                if J_type[j] == 'R':
                    JJ_t[i, j, :] = cross((A_I_j @ Ez), (RR[i, :] - RR[j, :] - A_I_j @ cc[j, j, :]))
                else:
                    JJ_t[i, j, :] = A_I_j @ Ez

                if BB[j] == -1:
                    j = -1
                else:
                    j = BB[j]
    return JJ_t


