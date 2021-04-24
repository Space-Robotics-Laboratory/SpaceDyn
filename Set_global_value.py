import numpy as np
import Global_Value as gl
import math





gl._init()

gl.set_value('num_q', 6)
gl.set_value('Gravity', np.array([0, 0, 0]))
gl.set_value('m0', 100)
gl.set_value('Ez', np.array([0, 0, 1]))
gl.set_value('m', np.array([10, 10, 10, 10, 10, 10]))
gl.set_value('inertia0', np.array([[10, 0, 0],
                                   [0, 10, 0],
                                   [0, 0, 10]]))

gl.set_value('inertia', np.array([[[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 0.1]],
                                  [[1, 0, 0],
                                   [0, 0.1, 0],
                                   [0, 0, 1]],
                                  [[1, 0, 0],
                                   [0, 0.1, 0],
                                   [0, 0, 1]],
                                  [[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 0.1]],
                                  [[1, 0, 0],
                                   [0, 0.1, 0],
                                   [0, 0, 1]],
                                  [[1, 0, 0],
                                   [0, 0.1, 0],
                                   [0, 0, 1]]]))

gl.set_value('d_time', 0.01)

gl.set_value('SS', np.array([[-1, 1, 0, 0, 0, 0],
                             [0, -1, 1, 0, 0, 0],
                             [0, 0, -1, 0, 0, 0],
                             [0, 0, 0, -1, 1, 0],
                             [0, 0, 0, 0, -1, 1],
                             [0, 0, 0, 0, 0, -1]]))

gl.set_value('S0', [1, 0, 0, 1, 0, 0])

gl.set_value('SE', [0, 0, 1, 0, 0, 1])

# gl.set_value('BB', [0, 1, 2, 0, 4, 5])

gl.set_value('BB', [-1, 0, 1, -1, 3, 4])

# gl.set_value('J_type', ['R', 'R', 'P', 'R', 'R', 'P'])
gl.set_value('J_type', ['R', 'R', 'R', 'R', 'R', 'R'])

gl.set_value('Qi', np.array([[-math.pi/2, 0, 0],
                             [math.pi/2, 0, 0],
                             [0, 0, 0],
                             [math.pi/2, 0, 0],
                             [math.pi/2, 0, 0],
                             [0, 0, 0]]))

gl.set_value('c0', np.array([[0, 1, 0], [0, 0, 0], [0, 0, 0], [0, -1, 0], [0, 0, 0], [0, 0, 0]]))

gl.set_value('ce', np.array([[0, 0, 0],
                             [0, 0, 0],
                             [0, 0.5, 0],
                             [0, 0, 0],
                             [0, 0, 0],
                             [0, 0.5, 0]]))

gl.set_value('Qe', np.array([[0, 0, 0],
                             [0, 0, 0],
                             [0, 0, math.pi/2],
                             [0, 0, 0],
                             [0, 0, 0],
                             [0, 0, math.pi/2]]))

gl.set_value('cc', np.zeros((6, 6, 3)))






'''

gl._init()

gl.set_value('num_q', 1)
gl.set_value('Gravity', np.array([0, 0, 0]))
gl.set_value('m0', 100)
gl.set_value('Ez', np.array([0, 0, 1]))
gl.set_value('m', np.array([100]))
gl.set_value('inertia0', np.array([[10, 0, 0],
                                   [0, 10, 0],
                                   [0, 0, 10]]))

gl.set_value('inertia', np.array([[[10, 0, 0],
                                   [0, 10, 0],
                                   [0, 0, 10]]]))

gl.set_value('mass', 200)
gl.set_value('d_time', 0.01)
gl.set_value('SS', np.array([[-1]]))

gl.set_value('S0', np.array([1]))

gl.set_value('SE', np.array([1]))

gl.set_value('BB', np.array([-1]))

gl.set_value('J_type', ['R'])

gl.set_value('Qi', np.array([[0, 0, 0]]))

gl.set_value('c0', np.array([[0, 0.5, 0]]))

gl.set_value('ce', np.array([[0, 0, 0.5]]))

gl.set_value('Qe', np.array([[0, 0, 0]]))

gl.set_value('cc', np.array([[[0, 0, -0.5]]]))

'''





'''
gl._init()

gl.set_value('num_q', 6)
gl.set_value('Gravity', np.array([0, 0, 0]))
gl.set_value('m0', 100)
gl.set_value('Ez', np.array([0, 0, 1]))
gl.set_value('m', np.array([10, 10, 10, 10, 10, 10]))
gl.set_value('inertia0', np.array([[10, 0, 0],
                                   [0, 10, 0],
                                   [0, 0, 10]]))

gl.set_value('inertia', np.array([[[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 0.1]],
                                  [[1, 0, 0],
                                   [0, 0.1, 0],
                                   [0, 0, 1]],
                                  [[1, 0, 0],
                                   [0, 0.1, 0],
                                   [0, 0, 1]],
                                  [[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 0.1]],
                                  [[1, 0, 0],
                                   [0, 0.1, 0],
                                   [0, 0, 1]],
                                  [[1, 0, 0],
                                   [0, 0.1, 0],
                                   [0, 0, 1]]]))

gl.set_value('d_time', 0.01)

gl.set_value('SS', np.array([[-1, 1, 0, 0, 0, 0],
                             [0, -1, 1, 0, 0, 0],
                             [0, 0, -1, 0, 0, 0],
                             [0, 0, 0, -1, 1, 0],
                             [0, 0, 0, 0, -1, 1],
                             [0, 0, 0, 0, 0, -1]]))

gl.set_value('S0', [1, 0, 0, 0, 0, 0])

gl.set_value('SE', [0, 0, 0, 0, 0, 1])

# gl.set_value('BB', [0, 1, 2, 0, 4, 5])

gl.set_value('BB', [-1, 0, 1, 2, 3, 4])

# gl.set_value('J_type', ['R', 'R', 'P', 'R', 'R', 'P'])
gl.set_value('J_type', ['R', 'R', 'R', 'R', 'R', 'R'])

gl.set_value('Qi', np.array([[-math.pi/2, 0, 0],
                             [math.pi/2, 0, 0],
                             [0, 0, 0],
                             [math.pi/2, 0, 0],
                             [math.pi/2, 0, 0],
                             [0, 0, 0]]))

gl.set_value('c0', np.array([[0, 1, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]))

gl.set_value('ce', np.array([[0, 0, 0],
                             [0, 0, 0],
                             [0, 0, 0],
                             [0, 0, 0],
                             [0, 0, 0],
                             [0, 0.5, 0]]))

gl.set_value('Qe', np.array([[0, 0, 0],
                             [0, 0, 0],
                             [0, 0, 0],
                             [0, 0, 0],
                             [0, 0, 0],
                             [0, 0, math.pi/2]]))

gl.set_value('cc', np.zeros((6, 6, 3)))





'''








