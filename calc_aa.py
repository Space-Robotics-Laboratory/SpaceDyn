import numpy as np
from Get_global_value import num_q
from Get_global_value import BB
from Get_global_value import J_type
from Get_global_value import Qi
from rpy2dc import rpy2dc


def calc_aa(A0, q):
    AA = np.zeros((num_q, 3, 3))
    if num_q == 0:
        print('Single body, there is no link')

    else:
        A_I_0 = A0

        for i in range(num_q):
            if BB[i] == -1:
                if J_type[i] == 'R':
                    A_0_i = rpy2dc(Qi[i, 0], Qi[i, 1], Qi[i, 2] + q[i]).T
                else:
                    A_0_i = rpy2dc(Qi[i, :], 0, 0).T

                AA[i, 0:3, 0:3] = A_I_0 @ A_0_i

            else:
                if J_type[i] == 'R':
                    A_BB_i = rpy2dc(Qi[i, 0], Qi[i, 1], Qi[i, 2] + q[i]).T
                else:
                    A_BB_i = rpy2dc(Qi[i, 0:3], 0, 0).T

                AA[i, 0:3, 0:3] = AA[BB[i], 0:3, 0:3] @ A_BB_i

    return AA
































