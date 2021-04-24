import numpy as np
from Cxyz import cxyz


def eul2dc(phi, theta, psi):
    C = cxyz(psi, 0, 0, 1) @ cxyz(theta, 1, 0, 0) @ cxyz(phi, 0, 0, 1)

    # if np.size(phi == 3):
    #     C = cxyz(phi[2], 0, 0, 1) @ cxyz(phi[1], 1, 0, 0) @ cxyz(phi[0], 0, 0, 1)
    # else:
    #     C = cxyz(psi, 0, 0, 1) @ cxyz(theta, 1, 0, 0) @ cxyz(phi, 0, 0, 1)

    return C
