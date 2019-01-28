# Import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import math

# Reading the input file name
br = 60 * "="
print(br + '\n' + 'PORTICO 2D' + '\n' + br + '\n' + '\n' + br + '\n' + br)
inputFileName = 'portico1'
fileName = inputFileName + '.fem'
with open(fileName, 'r') as myFile:
    file = myFile.read().splitlines()
myFile.close()


# Reading function
def readnum(keyword, file, type='int'):
    if keyword in file:
        totalnum = int(file[file.index(keyword) + 1])
        if totalnum > 0:
            num = pd.DataFrame(list(map(lambda x: x.split(),
                                        file[file.index(keyword) + 2:
                                             file.index(keyword) + 2 + totalnum])
                                    )
                               ).astype(type).as_matrix()
        else:
            num = pd.DataFrame()
        return totalnum, num
    else:
        raise ValueError("{} Not Found.".format(keyword))


def removefirstcol(df):
    return df[:, 1:]


# Reading nodal coordinates
print("READING NODAL COORDINATES...")
totalNumNodes, coords = readnum('*COORDINATES', file, type='float')
coords = removefirstcol(coords)

# Number of nodal degrees of freedom
numNodalDOFs = 3

# Reading element groups
print("READING ELEMENT GROUPS...")
numGroups, groups = readnum('*ELEMENT_GROUPS', file)
groups = removefirstcol(groups)

# Total number of elements
totalNumElements = groups.sum()

# Reading incidences
print("READING INCIDENCES...")
if '*INCIDENCES' in file:
    incid = pd.DataFrame(list(map(lambda groups: groups.split(),
                                  file[file.index('*INCIDENCES') + 1:
                                       file.index('*INCIDENCES') + totalNumElements + 1])
                              )
                         ).astype('int').as_matrix()
else:
    raise ValueError('{} Not Found.'.format('*INCIDENCES'))

# deletes element numbers
incid = removefirstcol(incid)

# Reading Materials
print('READING MATERIALS...')
numMaters, maters = readnum('*MATERIALS', file, type='float')

# Reading Geometric Properties
print('READING GEOMETRIC PROPERTIES...')
numGPs, gps = readnum('*GEOMETRIC_PROPERTIES', file, type='float')

# Reading Boundary Conditions
print('READING BOUNDARY CONDITIONS...')
numBCNodes, hdbcNodes = readnum('*BCNODES', file)

# Reading Nodal Loads
print('READING LOADS...')
numLoadedNodes, loads = readnum('*LOADS', file)
print(br)

# Reading Distributed Loads
print('READING DISTRIBUTED LOADS...')
numLoadedElements, distLoads = readnum('*DISTRIBUTED_LOADS', file)

# Element Orientation: length, lx = sin(theta), ly = cos(theta)
print('SOLUTION...\n' + br + '\n' + br)

# Element Lengths
lengths = np.sqrt(np.absolute((coords[incid[:, 0] - 1, 0] - coords[incid[:, 1] - 1, 0])**2 +
                              (coords[incid[:, 0] - 1, 1] - coords[incid[:, 1] - 1, 1])**2))

# lx = sin(theta)
lx = (coords[incid[:, 1] - 1, 1] - coords[incid[:, 0] - 1, 1]) / lengths

# ly = cos(theta)
ly = (coords[incid[:, 1] - 1, 0] - coords[incid[:, 0] - 1, 0]) / lengths

prop = np.matrix([lengths, lx, ly]).transpose()

# Distributed Load Global Matrix
globalDistLoads = np.zeros((totalNumElements, 2))

# Organize distributed loads in a matrix with all elements
globalDistLoads[distLoads[:, 0] - 1, :] = distLoads[:, 1:3]

# Solver
# DOF Numbering
# DOFs matriz initialized with ones
nodalDOFNumbers = np.ones((numNodalDOFs, totalNumNodes), dtype=int)

# Boundary Conditions on DOFs with Zero Values
for i in range(0, numBCNodes):
    nodalDOFNumbers[hdbcNodes[i, 1] - 1, hdbcNodes[i, 0] - 1] = 0

# Total Number of DOFs for the Model
totalNumDOFs = totalNumNodes * numNodalDOFs

# Numbers of Free and Restricted DOFs
totalNumFreeDOFs = np.sum(nodalDOFNumbers == 1)
totalNumRestrDOFs = np.sum(nodalDOFNumbers == 0)

# DOFs Numbering
ukeqnum = 0
keqnum = totalNumFreeDOFs
for i in range(totalNumNodes):
    for j in range(numNodalDOFs):
        if nodalDOFNumbers[j, i] == 1:
            ukeqnum = ukeqnum + 1
            nodalDOFNumbers[j, i] = ukeqnum
        elif nodalDOFNumbers[j, i] == 0:
            keqnum = keqnum + 1
            nodalDOFNumbers[j, i] = keqnum

# Assembling of the Global Stiffness Matrix and Global Distributed Loads Vector
# Creates Zero Global Matrix
kg = np.zeros((totalNumDOFs, totalNumDOFs))

# Distributed Loads
fgDistLoads = np.zeros((totalNumDOFs, 1))

elemNum = 0
for grp in range(0, numGroups):
    E = maters[grp, 0]
    rho = maters[grp, 3]
    A = gps[grp, 0]
    Iz = gps[grp, 1]

    for tne in range(0, groups[grp, 0]):
        c = prop[elemNum, 2]
        s = prop[elemNum, 1]
        l = prop[elemNum, 0]

        # Rotation Matrix
        T = np.array(([ c, s, 0,  0, 0, 0],
                      [-s, c, 0,  0, 0, 0],
                      [ 0, 0, 1,  0, 0, 0],
                      [ 0, 0, 0,  c, s, 0],
                      [ 0, 0, 0, -s, c, 0],
                      [ 0, 0, 0,  0, 0, 1]))

        # Element Stiffness Matrix
        # Axial Forces
        kb = E * A / l
        # Transversal Forces
        kf = E * Iz / l ** 3
        ke = np.array(([ kb,          0,               0, -kb,           0,               0],
                       [  0,    12 * kf,      6 * l * kf,   0,    -12 * kf,      6 * l * kf],
                       [  0, 6 * l * kf, 4 * l ** 2 * kf,   0, -6 * l * kf, 2 * l ** 2 * kf],
                       [-kb,          0,               0,  kb,           0,               0],
                       [  0,   -12 * kf,     -6 * l * kf,   0,     12 * kf,      -6 * l *kf],
                       [  0, 6 * l * kf, 2 * l ** 2 * kf,   0, -6 * l * kf, 4 * l ** 2 * kf]))
        ke = np.dot(np.dot((T.transpose()), ke), T)

        # Element DOFs
        elemEqs = nodalDOFNumbers[:, incid[elemNum, :] - 1]
        elemEqs = elemEqs.flatten('F').reshape(elemEqs.shape[0] * elemEqs.shape[1], 1)

        # Assembling into the Global Matrix
        kg[elemEqs - 1, elemEqs.reshape((-1)) - 1] = kg[elemEqs - 1, elemEqs.reshape((-1)) - 1] + ke

        # Distributed Loads
        # Applied Distributed Load Rotation Matrix
        Tq = np.array(([ c, s],
                       [-s, c]))

        # Element Distributed Load Vector
        q = np.dot(Tq, np.array((globalDistLoads[elemNum, 0],
                                 globalDistLoads[elemNum, 1] - rho * A * 9.81)))

        # Element Nodal Distributed Load Vector
        fe = l / 2 * np.array(([q[0], q[1], q[1] * l / 6, q[0], q[1], -q[1] * l / 6]))
        fe = np.dot(T.transpose(), fe)

        # Element DOFs
        elemEqs = nodalDOFNumbers[:, incid[elemNum, :] - 1]
        elemEqs = elemEqs.flatten('F')

        # Assembling into the Global Matrix
        fgDistLoads[elemEqs - 1, 0] = fgDistLoads[elemEqs - 1, 0] + fe

        elemNum += 1

# Assembling Global Load Vector
fgnLoads = np.zeros((totalNumDOFs, 1))

# Assembling Nodal Loads
for i in range(0, numLoadedNodes):
    # Element DOFs
    elemEqs = nodalDOFNumbers[loads[i, 1] - 1, loads[i, 0] - 1]

    # Assebmbling into the Global Vector
    fgnLoads[elemEqs - 1] = loads[i, 2]

fg = fgnLoads + fgDistLoads

# Solving system of equations
u = np.zeros((totalNumDOFs, 1))

u[0:totalNumFreeDOFs, :] = np.dot(np.linalg.inv(kg[:totalNumFreeDOFs,
                                                :totalNumFreeDOFs]),
                                  fg[:totalNumFreeDOFs, 0].reshape(totalNumFreeDOFs, 1))

# Calculates Reaction Forces, Element Internal Loads
# Reaction Forces
rf = (np.dot(kg[totalNumFreeDOFs:, 0:totalNumFreeDOFs], u[0:totalNumFreeDOFs]) -
      fg[totalNumFreeDOFs:])

# Deformed Coordinates
coordsd = coords + u[np.concatenate(nodalDOFNumbers[0:2, :]) - 1].reshape(coords.shape[1],
                                                                          coords.shape[0]).transpose()
scaleFactor = 8e2
coordsdd = coords + scaleFactor * u[np.concatenate(nodalDOFNumbers[0:2, :]) - 1].reshape(coords.shape[1],
                                                                                         coords.shape[0]).transpose()

# Element Normal Force, Shear Force, Internal Moment and Stresses for y-max and for y-min
nx = np.zeros((totalNumElements, 1))
vy = np.zeros((totalNumElements, 1))
mz = np.zeros((totalNumElements, 2))
stress_ymax = np.zeros((totalNumElements, 2))
stress_ymin = np.zeros((totalNumElements, 2))

# Calculates Strain and Stress for the Elements
elemNum = 0
for grp in range(0, numGroups):
    E = maters[grp, 0]
    A = gps[grp, 0]
    Iz = gps[grp, 1]
    ymax = gps[grp, 2]
    ymin = gps[grp, 3]

    for tne in range(0, groups[grp, 0]):
        c = prop[elemNum, 2]
        s = prop[elemNum, 1]
        l = prop[elemNum, 0]

        # Rotation Matrix
        T = np.array(([ c, s, 0,  0, 0, 0],
                      [-s, c, 0,  0, 0, 0],
                      [ 0, 0, 1,  0, 0, 0],
                      [ 0, 0, 0,  c, s, 0],
                      [ 0, 0, 0, -s, c, 0],
                      [ 0, 0, 0,  0, 0, 1]))

        # Element DOFs and Displacements
        elemEqs = nodalDOFNumbers[:, incid[elemNum, :] - 1]
        elemU = u[elemEqs[:] - 1].flatten('F')
        elemU = np.dot(T, elemU)

        # Normal Force
        nx[elemNum] = E * A * (-1 / l * elemU[0] + 1 / l * elemU[3])

        # Shear Force
        vy[elemNum] = E * Iz * 6 / l ** 2 * (2 / l * elemU[1] + elemU[2] -
                                             2 / l * elemU[4] + elemU[5])

        # Moments in Nodes 1 and 2
        mz[elemNum, 0] = E * Iz * (-6 / l ** 2 * elemU[1] - 4 / l * elemU[2] +
                                    6 / l ** 2 * elemU[4] - 2 / l * elemU[5])
        mz[elemNum, 1] = E * Iz * ( 6 / l ** 2 * elemU[1] + 2 / l * elemU[2] +
                                    -6 / l ** 2 * elemU[4] + 4 / l * elemU[5])

        # Stress
        stress_ymax[elemNum, 0] = nx[elemNum] / A - mz[elemNum, 0] / Iz * ymax
        stress_ymax[elemNum, 1] = nx[elemNum] / A - mz[elemNum, 1] / Iz * ymax
        stress_ymin[elemNum, 0] = nx[elemNum] / A - mz[elemNum, 0] / Iz * ymin
        stress_ymin[elemNum, 1] = nx[elemNum] / A - mz[elemNum, 1] / Iz * ymin

        elemNum += 1

# Output File
outFileName = inputFileName + 'py.out2'
with open(outFileName, 'w') as f:
    f.write('*DISPLACEMENTS\n')
    for i in range(totalNumNodes):
        f.write('{:d} {:.4f} {:.4f} {:.4f}\n'.format(i + 1, u[nodalDOFNumbers[0, i] - 1, 0],
                                                     u[nodalDOFNumbers[1, i] - 1, 0], u[nodalDOFNumbers[2, i] - 1, 0]))

    f.write('\n*REACTION_FORCES\n')
    for i in range(totalNumRestrDOFs):
        if hdbcNodes[i, 1] == 1:
            f.write('{:d} FX = {:.6e}\n'.format(hdbcNodes[i, 0], rf[i, 0]))
        elif hdbcNodes[i, 1] == 2:
            f.write('{:d} FY = {:.6e}\n'.format(hdbcNodes[i, 0], rf[i, 0]))
        else:
            f.write('{:d} MZ = {:.6e}\n'.format(hdbcNodes[i, 0], rf[i, 0]))

    f.write('\n*ELEMENT_NORMAL_FORCE\n')
    for i in range(totalNumElements):
        f.write('{:d} {:.6e}\n'.format(i + 1, nx[i].item()))

    f.write('\n*ELEMENT_SHEAR_FORCE\n')
    for i in range(totalNumElements):
        f.write('{:d} {:.6e}\n'.format(i + 1, vy[i].item()))

    f.write('\n*ELEMENT_MOMENT_LIMITS\n')
    for i in range(totalNumElements):
        f.write('{:d} {:.6e} {:.6e}\n'.format(i + 1, mz[i, 0], mz[i, 1]))

    f.write('\n*ELEMENT_STRESSES_YMAX\n')
    for i in range(totalNumElements):
        f.write('{:d} {:.6e} {:.6e}\n'.format(i + 1, stress_ymax[i, 0], stress_ymax[i, 1]))

    f.write('\n*ELEMENT_STRESSES_YMIN\n')
    for i in range(totalNumElements):
        f.write('{:d} {:.6e} {:.6e}\n'.format(i + 1, stress_ymin[i, 0], stress_ymin[i, 1]))

    f.close()

# Post Processing
print('POST-PROCESSING...\n' + br)


# Original Mesh Function
def originalMesh(ax, num_elements=totalNumElements):
    for elemNum in range(num_elements):
        ax.plot(coords[incid[elemNum, :] - 1, 0], coords[incid[elemNum, :] - 1, 1],
                color='white', linewidth=1, alpha=0.5)


# Nodes of the Underformed Mesh Function
def nodes_mesh(ax, coords=coords, order=1):
    ax.scatter(x=coords[:, 0], y=coords[:, 1], c='yellow',
               s=10, alpha=0.5, zorder=order)


# Circular Quiver Function
def plot_circ_quiver(h, k, r, M, color):
    if M > 0:
        t = np.linspace(-math.pi/3, 4 * math.pi/3)
    else:
        t = np.linspace(4 * math.pi/3, -math.pi/3)

    x = r * np.cos(t) + h
    y = r * np.sin(t) + k
    plt.plot(x, y, color=color)

    ax = plt.gca()

    if M < 0:
        ax.annotate('', xy=(x[-1] - r/3, y[-1]), xytext=(x[-2] - r/3, y[-2]),
                    arrowprops=dict(color=color, arrowstyle="->"))
    else:
        ax.annotate('', xy=(x[-1] + r/3, y[-1]), xytext=(x[-2] + r/3, y[-2]),
                    arrowprops=dict(color=color, arrowstyle="->"))


# Color gradient function
def plot_gradient_hack(p0, p1, npts=20, cmap=None, **kw):
    """
    Draw a gradient between p0 and p1 using a colormap
    The **kw dictionary gets passed to plt.plot, so things like linestyle,
    linewidth, labels, etc can be modified directly.
    """
    x_1, y_1 = p0
    x_2, y_2 = p1

    X = np.linspace(x_1, x_2, npts)
    Xs = X[:-1]
    Xf = X[1:]
    Xpairs = zip(Xs, Xf)

    Y = np.linspace(y_1, y_2, npts)
    Ys = Y[:-1]
    Yf = Y[1:]
    Ypairs = zip(Ys, Yf)

    C = np.linspace(0, 1, npts)
    cmap = plt.get_cmap(cmap)
    # the simplest way of doing this is to just do the following:
    for x, y, c in zip(Xpairs, Ypairs, C):
        plt.plot(x, y, '-', c=cmap(c), **kw)


def plot_gradient_rbg_pairs(p0, p1, rgb0, rgb1, **kw):
    """Form the gradient from RGB values at each point
    The **kw dictionary gets passed to plt.plot, so things like linestyle,
    linewidth, labels, etc can be modified directly.
    """
    cmap = LinearSegmentedColormap.from_list('jet', (rgb0, rgb1))
    plot_gradient_hack(p0, p1, cmap=cmap, **kw)


# Color definition Function
def rankmin(x):
    x = np.round(x, decimals=4)
    u, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    csum = np.zeros_like(counts)
    csum[1:] = counts[:-1].cumsum()
    return csum[inv]


# Mesh and Boundary Conditions
plt.style.use('dark_background')
fig1 = plt.figure()
plt.subplots_adjust(hspace=0.7, wspace=0.7)

ax1 = fig1.add_subplot(211)
ax1.set_title('Portico', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])

# Plots Original Mesh
originalMesh(ax1)

# Plots Nodes of the Underformed Mesh
nodes_mesh(ax1)

# Plot Loads
r = np.mean(prop[:, 0])
fFe = 0.25 * loads[:, 2] * abs(max(prop[:, 0]).item()) / abs(max(loads[:, 2]).item())
for i in range(loads.shape[0]):
    if loads[i, 1] == 1:
        ax1.quiver(coords[loads[i, 0] - 1, 0] - fFe[i], coords[loads[i, 0] - 1, 1], fFe[i], 0,
                   color='red', scale=1, units='xy', scale_units='xy', headlength=10,
                   headwidth=5)
    elif loads[i, 1] == 2:
        ax1.quiver(coords[loads[i, 0] - 1, 0], coords[loads[i, 0] - 1, 1] - fFe[i], 0, fFe[i],
                   color='red', scale=1, units='xy', scale_units='xy', headlength=10,
                   headwidth=5)
    elif loads[i, 1] == 3:
        plot_circ_quiver(coords[loads[i, 0] - 1, 0], coords[loads[i, 0] - 1, 1],
                         r/5, loads[i, 2], 'red')

# Plot Distributed Loads
distL = (0.15 * (distLoads[:, 1:3]) * np.asscalar(abs(max(prop[:, 0]))) /
         max(np.sqrt(distLoads[:, 1]**2 + distLoads[:, 2]**2)))
num_arrows_per_element = 11
for i in range(numLoadedElements):
    dx = ((coords[incid[distLoads[i, 0] - 1, 1] - 1, 0] -
          coords[incid[distLoads[i, 0] - 1, 0] - 1, 0])
          / (num_arrows_per_element - 1))
    dy = ((coords[incid[distLoads[i, 0] - 1, 1] - 1, 1] -
          coords[incid[distLoads[i, 0] - 1, 0] - 1, 1])
          / (num_arrows_per_element - 1))
    for j in range(num_arrows_per_element):
        ax1.quiver(coords[incid[distLoads[i, 0] - 1, 0] - 1, 0] - distL[i, 0] + j * dx,
                   coords[incid[distLoads[i, 0] - 1, 0] - 1, 1] - distL[i, 1] + j * dy, distL[i, 0], distL[i, 1],
                   color='red', scale=1, units='xy', scale_units='xy', headlength=10,
                   headwidth=5)

# Plot BCs
apoio = np.asscalar(0.7e-1*abs(max(prop[:, 0])))
for i in range(hdbcNodes.shape[0]):
    if hdbcNodes[i, 1] == 1:
        ax1.plot(coords[hdbcNodes[i, 0] - 1, 0] - apoio,
                 coords[hdbcNodes[i, 0] - 1, 1], marker='>', color='black',
                 markersize=8, markeredgecolor='blue')
    elif hdbcNodes[i, 1] == 2:
        ax1.plot(coords[hdbcNodes[i, 0] - 1, 0],
                 coords[hdbcNodes[i, 0] - 1, 1] - apoio,
                 marker='^', color='black', markersize=6, markeredgecolor='blue')
    elif hdbcNodes[i, 1] == 3:
        ax1.plot(coords[hdbcNodes[i, 0] - 1, 0], coords[hdbcNodes[i, 1] - 1, 1],
                 marker='X', color='blue', markersize=5)

# Plot Reaction Forces
ax2 = fig1.add_subplot(212)
ax2.set_title('Reaction Loads', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax2)

# Plot Reaction Loads
fRF = 0.25 * (rf * abs(max(prop[:, 0]))) / abs(max(rf[:, 0]))
for i in range(hdbcNodes.shape[0]):
    if hdbcNodes[i, 1] == 1:
        ax2.quiver(coords[hdbcNodes[i, 0] - 1, 0] - fRF[i, 0],
                   coords[hdbcNodes[i, 0] - 1, 1], fRF[i, 0], 0,
                   color='blue', scale=1, units='xy',
                   scale_units='xy', headlength=10, headwidth=5)
    elif hdbcNodes[i, 1] == 2:
        ax2.quiver(coords[hdbcNodes[i, 0] - 1, 0],
                   coords[hdbcNodes[i, 0] - 1, 1] - fRF[i, 0],
                   0, fRF[i, 0], color='blue', scale=1, units='xy',
                   scale_units='xy', headlength=10, headwidth=5)
    elif hdbcNodes[i, 1] == 3:
        plot_circ_quiver(coords[hdbcNodes[i, 0] - 1, 0], coords[hdbcNodes[i, 0] - 1, 1],
                         r / 5, fRF[i, 0], 'blue')

nodes_mesh(ax2, order=totalNumElements + 1)
fig1.show()

# Plot Element Stresses Ymin
fig2 = plt.figure()
plt.subplots_adjust(hspace=0.7, wspace=0.7)
ax3 = fig2.add_subplot(211)
ax3.set_title(r'$Element\ stresses\ y_{min}$', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax3)

# Plots Stresses on Elements of the Deformed Mesh
z = rankmin(stress_ymin[:, 0] - stress_ymin[:, 1])
colors = pl.cm.jet(np.linspace(0, 1, len(set(z))))
Z = np.unique(z, return_index=False)
for elemNum in range(totalNumElements):
    xy = coordsdd[incid[elemNum, :] - 1, :]
    plot_gradient_rbg_pairs(xy[0, :], xy[1, :], colors[np.where(z[incid[elemNum, 0] - 1] == Z)[0][0]],
                            colors[np.where(z[incid[elemNum, 1] - 1] == Z)[0][0]])

sm = plt.cm.ScalarMappable(cmap=plt.cm.jet,
                           norm=plt.Normalize(vmin=stress_ymin.min(), vmax=stress_ymin.max()))
sm._A = []
plt.colorbar(sm, ax=ax3, ticks=np.linspace(stress_ymin.min(), stress_ymin.max(), num=10))

nodes_mesh(ax3, coordsdd, order=totalNumElements + 1)

# Plot Element Stresses Ymax
ax4 = fig2.add_subplot(212)
ax4.set_title(r'$Element\ stresses\ y_{max}$', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax4)

# Plots Stresses on Elements of the Deformed Mesh
z2 = rankmin(stress_ymax[:, 0] - stress_ymax[:, 1])
colors2 = pl.cm.jet(np.linspace(0, 1, len(set(z2))))
Z2 = np.unique(z2, return_index=False)
for elemNum in range(totalNumElements):
    xy = coordsdd[incid[elemNum, :] - 1, :]
    plot_gradient_rbg_pairs(xy[0, :], xy[1, :], colors2[np.where(z2[incid[elemNum, 0] - 1] == Z)[0][0]],
                            colors2[np.where(z2[incid[elemNum, 1] - 1] == Z2)[0][0]])

sm = plt.cm.ScalarMappable(cmap=plt.cm.jet,
                           norm=plt.Normalize(vmin=stress_ymax.min(), vmax=stress_ymax.max()))
sm._A = []
plt.colorbar(sm, ax=ax4, ticks=np.linspace(stress_ymax.min(), stress_ymax.max(), num=10))

nodes_mesh(ax4, coordsdd, order=totalNumElements + 1)
fig2.show()

# Plot Element Normal Force
fig3 = plt.figure()
plt.subplots_adjust(hspace=0.7, wspace=0.7)
ax5 = fig3.add_subplot(211)
ax5.set_title(r'$Element\ normal\ forces\ N_{x}$', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax5)

# Plot Stresses on Elements of the Deformed Mesh
colors = plt.cm.jet(nx / float(nx.max()))
for elemNum in range(totalNumElements):
    ax5.plot(coordsdd[incid[elemNum, ] - 1, 0],
             coordsdd[incid[elemNum, ] - 1, 1],
             color=colors[elemNum].flatten())

sm = plt.cm.ScalarMappable(cmap=plt.cm.jet,
                           norm=plt.Normalize(vmin=nx.min(), vmax=nx.max()))
sm._A = []
plt.colorbar(sm, ax=ax5, ticks=np.linspace(nx.min(), nx.max(), num=10))
nodes_mesh(ax5, coordsdd, order=totalNumElements + 1)

# Plot Elements Shear Force
plt.subplots_adjust(hspace=0.7, wspace=0.7)
ax6 = fig3.add_subplot(212)
ax6.set_title(r'$Element\ shear\ forces\ V_{y}$', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax6)

# Plot Stresses on Elements of the Deformed Mesh
colors = plt.cm.jet(vy / float(vy.max()))
for elemNum in range(totalNumElements):
    ax6.plot(coordsdd[incid[elemNum, ] - 1, 0],
             coordsdd[incid[elemNum, ] - 1, 1],
             color=colors[elemNum].flatten())

sm = plt.cm.ScalarMappable(cmap=plt.cm.jet,
                           norm=plt.Normalize(vmin=vy.min(), vmax=vy.max()))

sm._A = []
plt.colorbar(sm, ax=ax6, ticks=np.linspace(vy.min(), vy.max(), num=10))
nodes_mesh(ax6, coordsdd, order=totalNumElements + 1)

plt.show()
