import numpy as np
from Get_global_value import d_time
from f_dyn import f_dyn
from skew_sym import skew_sym
from aw import aw


def f_dyn_rk2(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau):

    C0 = A0.T

    # 1st step
    tmp_vd0, tmp_wd0, tmp_qdd = f_dyn(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau)

    k1_R0 = d_time * v0
    k1_C0 = d_time * skew_sym(w0).T @ C0
    k1_A0 = k1_C0.T
    k1_q = d_time * qd

    k1_v0 = d_time * tmp_vd0
    k1_w0 = d_time * tmp_wd0
    k1_qd = d_time * tmp_qdd

    # 2nd step
    tmp_vd0, tmp_wd0, tmp_qdd = f_dyn(R0 + k1_R0/2, A0 + k1_A0/2, v0 + k1_v0/2, w0 + k1_w0/2, q + k1_q/2, qd + k1_qd/2, F0, T0, Fe, Te, tau)

    k2_R0 = d_time * (v0 + k1_v0 / 2)
    k2_C0 = d_time * (skew_sym(w0 + k1_w0 / 2).T @ C0)
    k2_A0 = k2_C0.T
    k2_q = d_time * (qd + k1_qd / 2)

    k2_v0 = d_time * tmp_vd0
    k2_w0 = d_time * tmp_wd0
    k2_qd = d_time * tmp_qdd

    # 3rd step
    tmp_vd0, tmp_wd0, tmp_qdd = f_dyn(R0 + k2_R0/2, A0 + k2_A0/2, v0 + k2_v0/2, w0 + k2_w0/2, q + k2_q/2, qd + k2_qd/2, F0, T0, Fe, Te, tau)

    k3_R0 = d_time * (v0 + k2_v0 / 2)
    k3_C0 = d_time * (skew_sym(w0 + k2_w0 / 2).T @ C0)
    k3_A0 = k3_C0.T
    k3_q = d_time * (qd + k2_qd / 2)

    k3_v0 = d_time * tmp_vd0
    k3_w0 = d_time * tmp_wd0
    k3_qd = d_time * tmp_qdd

    # 4th step
    tmp_vd0, tmp_wd0, tmp_qdd = f_dyn(R0 + k3_R0, A0 + k3_A0, v0 + k3_v0, w0 + k3_w0, q + k3_q, qd + k3_qd, F0, T0, Fe, Te, tau)
    k4_R0 = d_time * (v0 + k3_v0)
    k4_C0 = d_time * (skew_sym(w0 + k3_w0).T @ C0)
    k4_A0 = k4_C0.T
    k4_q = d_time * (qd + k3_qd)

    k4_v0 = d_time * tmp_vd0
    k4_w0 = d_time * tmp_wd0
    k4_qd = d_time * tmp_qdd

    R0_next = R0 + (k1_R0 + 2 * k2_R0 + 2 * k3_R0 + k4_R0) / 6
    C0_next = C0 + (k1_C0 + 2 * k2_C0 + 2 * k3_C0 + k4_C0) / 6
    A0_next = C0_next.T
    q_next = q + (k1_q + 2 * k2_q + 2 * k3_q + k4_q) / 6

    v0_next = v0 + (k1_v0 + 2 * k2_v0 + 2 * k3_v0 + k4_v0) / 6
    w0_next = w0 + (k1_w0 + 2 * k2_w0 + 2 * k3_w0 + k4_w0) / 6
    qd_next = qd + (k1_qd + 2 * k2_qd + 2 * k3_qd + k4_qd) / 6

    R0 = R0_next
    A0 = A0_next
    q = q_next
    v0 = v0_next
    w0 = w0_next
    qd = qd_next

    return R0, A0, v0, w0, q, qd
