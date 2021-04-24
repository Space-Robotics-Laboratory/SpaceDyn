import numpy as np
from math import cos
from math import sin


def cxyz(theta, x, y, z):
    direction_cosines = np.zeros((3, 3))
    if x == 1:
        direction_cosines = np.array([[1, 0, 0],
                                      [0, cos(theta), sin(theta)],
                                      [0, -sin(theta), cos(theta)]])
    elif y == 1:
        direction_cosines = np.array([[cos(theta), 0, -sin(theta)],
                                      [0, 1, 0],
                                      [sin(theta), 0, cos(theta)]])
    elif z == 1:
        direction_cosines = np.array([[cos(theta), sin(theta), 0],
                                      [-sin(theta), cos(theta), 0],
                                      [0, 0, 1]])

    return direction_cosines
