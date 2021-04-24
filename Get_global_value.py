import Global_Value as gl
import numpy as np


num_q = gl.get_value('num_q')
m0 = gl.get_value('m0')
m = gl.get_value('m')


mass = sum(m)+m0

# print(mass)

Ez = gl.get_value('Ez')
inertia = gl.get_value('inertia')
inertia0 = gl.get_value('inertia0')
BB = gl.get_value('BB')
J_type = gl.get_value('J_type')
Qi = gl.get_value('Qi')
Qe = gl.get_value('Qe')
c0 = gl.get_value('c0')

cc = gl.get_value('cc')

cc[0, 0, :] = np.array([0, 0, -0.5])
cc[1, 1, :] = np.array([0, 0, -0.5])
cc[2, 2, :] = np.array([0, 0, -0.5])
cc[3, 3, :] = np.array([0, 0, -0.5])
cc[4, 4, :] = np.array([0, 0, -0.5])
cc[5, 5, :] = np.array([0, 0, -0.5])
cc[0, 1, :] = np.array([0, 0, 0.5])
cc[1, 2, :] = np.array([0, 0, 0.5])
cc[3, 4, :] = np.array([0, 0, 0.5])
cc[4, 5, :] = np.array([0, 0, 0.5])


ce = gl.get_value('ce')
SS = gl.get_value('SS')
SE = gl.get_value('SE')
S0 = gl.get_value('S0')
Gravity = gl.get_value('Gravity')
d_time = gl.get_value('d_time')
