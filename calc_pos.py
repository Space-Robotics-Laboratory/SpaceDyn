import numpy as np
from Get_global_value import num_q
from Get_global_value import BB
from Get_global_value import J_type
from Get_global_value import Qi
from Get_global_value import c0
from Get_global_value import cc
from Get_global_value import Ez
from rpy2dc import rpy2dc


def calc_pos(R0, A0, AA, q):
    RR = np.zeros((num_q, 3))
    if num_q == 0:
        print('Single body, there is no link')
    else:
        for i in range(num_q):

            A_I_i = AA[i, 0:3, 0:3]

            if BB[i] == -1:
                if J_type[i] == 'R':
                    RR[i, 0:3] = R0[0:3] + A0 @ c0[i, 0:3] - A_I_i @ cc[i, i, 0:3]
                else:
                    RR[i, 0:3] = R0[0:3] + A0 @ c0[i, 0:3] + A_I_i @ (q[i] * Ez - cc[i, i, 0:3])

            else:
                A_I_BB = AA[BB[i], 0:3, 0:3]

                if J_type[i] == 'R':
                    RR[i, :] = RR[BB[i], :] + A_I_BB @ cc[BB[i], i, :] - A_I_i @ cc[i, i, :]
                else:
                    RR[i, :] = RR[BB[i], :] + A_I_BB @ cc[BB[i], i, :] + A_I_i @ (q[i] * Ez - cc[i, i, 0:3])

    return RR





















































