import numpy as np
from scipy.linalg import eigh
import math
from matplotlib import pyplot as plt
import femeq_1d as fem

# Input Parameters Data

# Length of the bar - metre:
lenx1 = float(input('Enter the length of the bar in m: '))
# Young's Modulus N/m2:
E = float(input('Enter the Modulous of Elasticity in N/m2 : '))
# Cross section in m2:
csareal1 = float(input('Enter the Area of cross section in m2: '))

a1 = float(input('Enter the first constant of linear f: '))
b1 = float(input('Enter the second constant of linear f: '))

a2 = float(input('Enter the first constant of quadratic f: '))
b2 = float(input('Enter the second constant of quadratic f: '))
c2 = float(input('Enter the third constant of quadratic f: '))

# Generate nodes, elements and connectivity

# Number of elements:
nelem1 = int(input('Enter the total number of Elements : '))
nnodes = nelem1 + 1  # number of nodes
coordx = np.linspace(0, lenx1, nnodes)

# Define the connectivity table:
connx = np.zeros((nelem1, 2))
for ne in range(0, nelem1, 1):
    connx[ne, :] = np.array([ne, ne+1])
# end of for loop
del ne

connxq = np.zeros((nelem1, 3))
for ne in range(0, nelem1 - 1, 1):
    connxq[ne, :] = np.array([ne, ne+1, ne+2])
# end of for loop
del ne

# Evaluate the stiffness matrix
matkg = fem.stiffmat(E, csareal1, nnodes, nelem1, connx, coordx)

matkgq = fem.stiffmatq(E, csareal1, nnodes, nelem1, connxq, coordx)

# Evaluate the body force vector
vecfg = fem.body_force(csareal1, nnodes, nelem1, connx, coordx, a1, b1)

vecfgq = fem.body_forceq(csareal1, nnodes, nelem1, connxq, coordx, a2, b2, c2)


# ----------------------Boundry condition and Loading---------------------#

displist = []
forcelist = []
for i in range(nnodes):
    a = str('u')+str(i+1)
    displist.append(a)
    c = str('fx')+str(i+1)
    forcelist.append(c)

forcemat = np.zeros((nnodes, 1))
tlon = int(input('Enter the total number of loaded nodes : ')) #total number of loaded nodes
for i in range(tlon):
    lon = int(input('\nEnter the node number of Loading : ')) #Loaded node
    fx = float(input('Enter the load at this node in N : '))
    forcemat[lon-1, 0] = fx
    vecfg = vecfg + forcemat
    vecfgq = vecfgq + forcemat

dispmat = np.ones((nnodes, 1))
tsupn = int(input('Enter the total number of nodes having supports : ')) #total number of supported nodes
for i in range(tsupn):
    supn = int(input('\nEnter the node number of support : ')) #supported node
    dispmat[supn-1, 0] = 0
    matkg[supn - 1, :] = 0
    matkg[supn - 1, supn - 1] = 1 / E
    matkgq[supn -1, :] = 0
    matkgq[supn - 1, supn - 1] = 1 / E
    vecfg[supn - 1, 0] = 0
    vecfgq[supn - 1, 0] = 0




# Solve the finite element equation
dispv = np.linalg.solve(matkg, vecfg)
dispvq = np.linalg.solve(matkgq, vecfgq)

# Deformation from the FEM Solution
ufe = np.zeros(nnodes)
for ne in range(0, nnodes, 1):
    ufe[ne] = dispv[ne]
# end of for loop
del ne

ufeq = np.zeros(nnodes)
for ne in range(0, nnodes, 1):
    ufeq[ne] = dispvq[ne]
# end of for loop
del ne


# Plotting the solutions

plt.plot(coordx, ufe, 'ro-')
plt.plot(coordx, ufeq, 'bo-')
plt.title('Deformation plot')
plt.xlabel('distance')
plt.ylabel('deformation')
plt.grid()
plt.show()

