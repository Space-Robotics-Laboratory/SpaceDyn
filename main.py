

import Global_Value
import Set_global_value
import Get_global_value

from Get_global_value import d_time
from Get_global_value import num_q
import numpy as np
import math
from eul2dc import eul2dc
from f_dyn_rk2 import f_dyn_rk2
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm




q = np.zeros(num_q)
qd = np.zeros(num_q)
qdd = np.zeros(num_q)

vv = np.zeros((num_q, 3))
ww = np.zeros((num_q, 3))
vd = np.zeros((num_q, 3))
wd = np.zeros((num_q, 3))

v0 = np.array([0, 0, 0])
w0 = np.array([0, 0, 0])
vd0 = np.array([0, 0, 0])
wd0 = np.array([0, 0, 0])

R0 = np.array([0, 0, 0])
Q0 = np.array([0, 0, 0])
A0 = np.eye(3)

Fe = np.zeros((num_q, 3))
Te = np.zeros((num_q, 3))
F0 = np.array([0, 0, 0])
T0 = np.array([0, 0, 0])

tau = np.zeros(num_q)

# print('ok1')

q1tempt = []
qdtempt = []
v0tempt = []
w0tempt = []
timetempt = []

# PID parameters set
desired_q = np.array([0.3, 0.2, 0.1, 0.6, 0.5, 0.4])
gain_spring = 10
gain_dumper = 10

for time in np.arange(0, 20, d_time):

    timetempt.append(time)

    # tau = np.array([0, 0, 0, 0, 0, 0])
    tau = gain_spring * (desired_q - q) - gain_dumper * qd

    R0, A0, v0, w0, q, qd = f_dyn_rk2(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau)

    # print("wait %f" % time)

    qdtempt.append(qd[0])
    q1tempt.append(q[3])
    v0tempt.append(v0[0])
    w0tempt.append(w0[0])

plt.plot(timetempt, q1tempt, linewidth=1.0, color='red', linestyle='--', label='q1tempt')
# plt.plot(timetempt, qdtempt, linewidth=1.0, color='black', linestyle='-.', label='qdtempt')
# plt.plot(timetempt, v0tempt, linewidth=1.0, color='green', linestyle='-', label='v0tempt')
# plt.plot(timetempt, w0tempt, linewidth=1.0, color='blue', linestyle=':', label='w0tempt')

# plt.legend(loc='upper left')

plt.grid(True)

plt.show()








'''
import test2
import test3
import test4
from test4 import ROOT

print(ROOT)

'''





'''
# use Monte Carlo to generate working space 
from calc_aa import calc_aa
from j_num import j_num
from calc_pos import calc_pos
from Get_global_value import SE
from f_kin_e import f_kin_e


A0 = np.eye(3)
R0 = np.array([0, 0, 0])

# 在[1, 10)之间均匀抽样，数组形状(1,6)
# q = np.random.uniform(-170/180 * math.pi, 170/180 * math.pi, (6,))
# print(q)
i = 0
POS_e1 = np.zeros((2000, 3))
ORI_e1 = np.zeros((2000, 3, 3))
POS_e2 = np.zeros((2000, 3))
ORI_e2 = np.zeros((2000, 3, 3))

while i < 2000:
    # A0 = np.random.rand(3, 3)

###############################################################################################################
    phi = np.random.uniform(-math.pi, math.pi)
    theta = np.random.uniform(-math.pi, math.pi)
    psi = np.random.uniform(-math.pi, math.pi)
    A0 = eul2dc(phi, theta, psi)
###############################################################################################################
    R0 = np.array([0, 0, 0])
    

    q = np.random.uniform(-170 / 180 * math.pi, 170 / 180 * math.pi, (6,))
    AA = calc_aa(A0, q)
    RR = calc_pos(R0, A0, AA, q)

    joints = j_num(0)
    POS_e1[i, :], ORI_e1[i, :, :] = f_kin_e(RR, AA, joints)

    # joints = j_num(1)
    # POS_e2[i, :], ORI_e2[i, :, :] = f_kin_e(RR, AA, joints)

    i += 1


fig = plt.figure()
ax = Axes3D(fig)

ax.scatter(POS_e1[:, 0], POS_e1[:, 1], POS_e1[:, 2])
# ax.scatter(POS_e2[:, 0], POS_e2[:, 1], POS_e2[:, 2])



# 添加坐标轴(顺序是Z, Y, X)
ax.set_zlabel('Z', fontdict={'size': 15, 'color': 'red'})
ax.set_ylabel('Y', fontdict={'size': 15, 'color': 'red'})
ax.set_xlabel('X', fontdict={'size': 15, 'color': 'red'})


##############################################################################################################
X = POS_e1[:, 0]
Y = POS_e1[:, 1]
Z = POS_e1[:, 2]
# ax.scatter(X, Y, Z, 'b-', linewidth=4, label='curve')

null = [6]*len(Z)
ax.scatter(null, Y, Z)
ax.scatter(X, null, Z)
ax.scatter(X, Y, null)
#############################################################################################################



plt.show()
'''





















