import numpy as np
from Cxyz import cxyz


def rpy2dc(roll, pitch, yaw):
    if np.size(roll) == 3:
        direction_cosine = cxyz(roll[2], 0, 0, 1) @ cxyz(roll[1], 0, 1, 0) @ cxyz(roll[0], 1, 0, 0)
    else:
        direction_cosine = cxyz(yaw, 0, 0, 1) @ cxyz(pitch, 0, 1, 0) @ cxyz(roll, 1, 0, 0)

    return direction_cosine
