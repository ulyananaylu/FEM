# Python module for 2-nodes bar elements
import numpy as np

# Stiffness Matrix Evaluation:

def stiffmat(Ex, areax, nnodes, nelem, conn, coord):
    matkg = np.zeros((nnodes, nnodes))  # global stiffness matrix
    for ne in range(0, nelem, 1):
        nn1 = int(conn[ne, 0])
        nn2 = int(conn[ne, 1])
        lenx = coord[nn2] - coord[nn1]
        conterm = areax * Ex / lenx
        # stiffness matrix for one element:
        matkl = conterm*np.array([[1, -1], [-1, 1]])
        # Assembly the global stiffness matrix & Body force
        nv = [nn1, nn2]
        nc = 0
        # Assembling the global stiffness matrix
        for jj in nv:
            matkg[jj, nv] = matkg[jj, nv] + matkl[nc, :]
            nc = nc + 1
    return matkg
# End of Stiffness matrix

# Force vector for body force
def body_force(areax, nnodes, nelem, conn, coord, a, b):
    vecfg = np.zeros((nnodes, 1))  # global force vector
    for ne in range(0, nelem, 1):
        nn1 = int(conn[ne, 0])
        nn2 = int(conn[ne, 1])
        lenx = coord[nn2] - coord[nn1]
        const1 = areax * lenx/2
        vecfb = const1 * np.array([[(a * coord[nn1] + b) * 1], [(a * coord[nn2] + b) * 1]])
        nv = [nn1, nn2]
        nc = 0
        for jj in nv:
            vecfg[jj, 0] = vecfg[jj, 0] + vecfb[nc, 0]
            nc = nc + 1
    return vecfg
# End of body force vector term


# Stiffness Quadratic Matrix Evaluation:

def stiffmatq(Ex, areax, nnodes, nelem, conn, coord):
    matkg = np.zeros((nnodes, nnodes))  # global stiffness matrix
    for ne in range(0, nelem-1, 1):
        nn1 = int(conn[ne, 0])
        nn2 = int(conn[ne, 1])
        nn3 = int(conn[ne, 2])
        lenx = coord[nn3] - coord[nn2]
        conterm = areax * Ex / lenx
        # stiffness matrix for one element:
        matkl = conterm*np.array([[7 / 6, -4 / 3, 1 / 6], [-4 / 3, 8 / 3, - 4 / 3], [1 / 6, -4 / 3, 7 / 6]])
        # Assembly the global stiffness matrix & Body force
        nv = [nn1, nn2, nn3]
        nc = 0
        # Assembling the global stiffness matrix
        for jj in nv:
            matkg[jj, nv] = matkg[jj, nv] + matkl[nc, :]
            nc = nc + 1
    return matkg
# End of Stiffness matrix


# Quadratic Force vector for body force
def body_forceq(areax, nnodes, nelem, conn, coord, a, b, c):
    vecfg = np.zeros((nnodes, 1))  # global force vector
    for ne in range(0, nelem-1, 1):
        nn1 = int(conn[ne, 0])
        nn2 = int(conn[ne, 1])
        nn3 = int(conn[ne, 2])
        lenx = coord[nn3] - coord[nn2]
        const1 = areax * lenx/2
        vecfb = const1 * np.array([[(a * (coord[nn1]**2) + b * coord[nn1] + c) * 1 / 3], [(a * (coord[nn2]**2) + b * coord[nn1] + c) * 4 / 3], [(a * (coord[nn3]**2) + b * coord[nn3] + c)* 1 / 3]])
        nv = [nn1, nn2, nn3]
        nc = 0
        for jj in nv:
            vecfg[jj, 0] = vecfg[jj, 0] + vecfb[nc, 0]
            nc = nc + 1
    return vecfg
# End of body force vector term
