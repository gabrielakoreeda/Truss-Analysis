
# Importing libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm

# Reading the input filename
br = 60 * '='
print(br + '\n' + 'TRUSS 2D' + '\n' + br + '\n' + '\n' + br + '\n' + br)
inpFileName = 'example2_truss'
fileName = inpFileName + '.fem'
with open(fileName, 'r') as myFile:
    file = myFile.read().splitlines()
myFile.close()


# Reading function
def readnum(keyword, document, type_num='int'):
    if keyword in document:
        totalnum = int(document[document.index(keyword) + 1])
        if totalnum > 0:
            num = pd.DataFrame(list(map(lambda x: x.split(),
                                        document[document.index(keyword) + 2:
                                                 document.index(keyword) + 2 + totalnum])
                                    )
                               ).astype(type_num).as_matrix()
        else:
            num = pd.DataFrame()
        return totalnum, num
    else:
        raise ValueError("{} Not Found.".format(keyword))


def removefirstcol(df):
    return df[:, 1:]


# Reading nodal coordinates
print('READING NODAL COORDINATES...')
totalNumNodes, coords = readnum('*COORDINATES', file, type_num='float')
coords = removefirstcol(coords)

# Number of nodal dofs
numNodalDOFs = 2

# Reading element groups
print('READING ELEMENT GROUPS...')
numGroups, groups = readnum('*ELEMENT_GROUPS', file)
groups = removefirstcol(groups)

# Total number of elements
totalNumElements = groups.sum()

# Reading incidences
print('READING ELEMENT INCIDENCES...')
if '*INCIDENCES' in file:
    incid = pd.DataFrame(list(map(lambda groups: groups.split(),
                                  file[file.index('*INCIDENCES') + 1:
                                       file.index('*INCIDENCES') + totalNumElements + 1])
                              )).astype('int').as_matrix()
else:
    raise ValueError('{} Not Found'.format('*INCIDENCES'))

# deletes element numbers
incid = removefirstcol(incid)

# Reading Materials
print('READING MATERIALS...')
numMaters, maters = readnum('*MATERIALS', file, type_num='float')

# Reading geometric properties
print('READING GEOMETRIC PROPERTIES...')
numGPs, gps = readnum('*GEOMETRIC_PROPERTIES', file, type_num='float')

# Reading boundary conditions
print('READING BOUNDARY CONDITIONS...')
numBCNodes, hdbcNodes = readnum('*BCNODES', file)

# Reading loads
print('READING LOADS...')
numLoadedNodes, loads = readnum('*LOADS', file)
print(br)

# Element orientation: length, lx = sin theta, ly = cos theta
print('SOLUTION...\n' + br + '\n' + br)

# Solver
# DOF Numbering
# DOFs matrix initialized with ones
nodalDOFNumbers = np.ones((numNodalDOFs, totalNumNodes), dtype=int)

# Boundary conditions on dofs with zero values
nodalDOFNumbers[hdbcNodes[:numBCNodes, 1] - 1, hdbcNodes[:numBCNodes, 0] - 1] = 0

# Total number of DOFs for the model
totalNumDOFs = totalNumNodes*numNodalDOFs

# Number of free and restricted dofs
totalNumFreeDOFs = np.sum(nodalDOFNumbers == 1)
totalNumRestrDOFs = np.sum(nodalDOFNumbers == 0)

# DOFs numbering
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

# element lengths
lengths = np.sqrt(np.absolute((coords[incid[:, 0] - 1, 0] - coords[incid[:, 1] - 1, 0])**2 +
                              (coords[incid[:, 0] - 1, 1] - coords[incid[:, 1] - 1, 1])**2))

# lx = sin
lx = (coords[incid[:, 1] - 1, 1] - coords[incid[:, 0] - 1, 1]) / lengths

# ly = cos
ly = (coords[incid[:, 1] - 1, 0] - coords[incid[:, 0] - 1, 0]) / lengths

prop = np.matrix([lengths, lx, ly]).transpose()

# Assembling of the global stiffness matrix and load vector
# Creates zero global matrix
kg = np.zeros((totalNumDOFs, totalNumDOFs))
elemNum = 0
for grp in range(numGroups):
    e = maters[grp, 0]
    a = gps[grp, 0]
    for tne in range(groups[grp - 1, 0]):
        cc = prop[elemNum, 2]**2
        ss = prop[elemNum, 1]**2
        cs = prop[elemNum, 1]*prop[elemNum, 2]

        # Element stiffness matrix
        kedf = [[cc, cs, -cc, -cs],
                [cs, ss, -cs, -ss],
                [-cc, -cs, cc, cs],
                [-cs, -ss, cs, ss]]
        ke = np.dot((e * a / prop[elemNum, 0]), kedf)

        # Element DOFs
        elemEqs = nodalDOFNumbers[:, incid[elemNum, :] - 1]

        # Assembling into the global matrix
        index = np.concatenate(np.column_stack(elemEqs)) - 1
        kg[index, index.reshape((-1, 1))] = kg[index, index.reshape((-1, 1))] + ke
        elemNum = elemNum + 1

# Assembling global matrix
fg = np.zeros((totalNumDOFs, 1))
for i in range(numLoadedNodes):
    # element dofs
    elemEqs = nodalDOFNumbers[loads[i, 1] - 1, loads[i, 0] - 1]

    # assembling into the global vector
    fg[elemEqs - 1, 0] = loads[i, 2]

# Solving systems of equations
u = np.zeros((totalNumDOFs, 1))
u[:totalNumFreeDOFs, 0] = np.dot(np.linalg.inv(kg[:totalNumFreeDOFs, :totalNumFreeDOFs]),
                                 fg[:totalNumFreeDOFs, 0])

# Calculates reaction forces, element strains and stresses
# reaction forces
rf = np.dot(kg[totalNumFreeDOFs:, :totalNumFreeDOFs], u[:totalNumFreeDOFs])

# Deformed coordinates
coordsd = coords + u[np.concatenate(np.row_stack(nodalDOFNumbers)) - 1].reshape(
    (nodalDOFNumbers.shape[0], nodalDOFNumbers.shape[1])).transpose()
scaleFactor = 100
coordsdd = coords + scaleFactor*u[np.concatenate(np.row_stack(nodalDOFNumbers)) - 1].reshape(
    (nodalDOFNumbers.shape[0], nodalDOFNumbers.shape[1])).transpose()

# Element strain and stress vectors
strain = np.zeros((totalNumElements, 1))
stress = np.zeros((totalNumElements, 1))

# Calculates strain and stress for the elements
elemNum = 0
for grp in range(numGroups):
    e = maters[grp, 0]
    for elNum in range(groups[grp - 1, 0]):
        # element dofs and displacements
        elemEqs = nodalDOFNumbers[:, incid[elemNum, :] - 1]
        index = np.concatenate(np.column_stack(elemEqs)) - 1
        elemU = u[index, 0]
        # strain
        strain[elemNum] = 1 / prop[elemNum, 0]*np.dot(
            np.array([[-prop[elemNum, 2], -prop[elemNum, 1],
                       prop[elemNum, 2], prop[elemNum, 1]]]), elemU)

        # stress
        stress[elemNum] = strain[elemNum]*e

        elemNum = elemNum + 1

# Post-Processing
print('POST-PROCESSING... \n' + br)

# Original Mesh Function
def originalMesh(ax, num_elements=totalNumElements):
    for elemNum in range(num_elements):
        ax.plot(coords[incid[elemNum, :] - 1, 0], coords[incid[elemNum, :] - 1, 1],
                color='white', linewidth=1, alpha=0.5)


# Nodes of the Underformed Mesh Function
def nodes_mesh(ax, coords=coords, order=1):
    ax.scatter(x=coords[:, 0], y=coords[:, 1], c='yellow',
               s=10, alpha=0.5, zorder=order)


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

# Mesh and boundary conditions
plt.style.use('dark_background')
fig1 = plt.figure()
plt.subplots_adjust(hspace=0.7, wspace=0.7)

ax1 = fig1.add_subplot(211)
ax1.set_title('Truss', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])

# Plots Original Mesh
originalMesh(ax1)

# Plots Nodes of the Undeformed Mesh
nodes_mesh(ax1)

# Plot loads
fFe = 0.25 * loads[:, 2] * np.asscalar(abs(max(prop[:, 0]))/abs(max(loads[:, 2])))
for i in range(len(loads)):
    if loads[i, 1] == 1:
        ax1.quiver(coords[loads[i, 0] - 1, 0] - fFe[i], coords[loads[i, 0] - 1, 1],
                   fFe[i], 0, color='red', scale=1, units='xy', scale_units='xy',
                   headlength=10, headwidth=5)
    elif loads[i, 1] == 2:
        ax1.quiver(coords[loads[i, 0] - 1, 0], coords[loads[i, 0] - 1, 1] - fFe[i],
                   0, fFe[i], color='red', scale=1, units='xy', scale_units='xy',
                   headlength=10, headwidth=5)

# Plot bcs
apoio = np.asscalar(0.7e-1*abs(max(prop[:, 0])))
for i in range(len(hdbcNodes)):
    if hdbcNodes[i, 1] == 1:
        ax1.plot(coords[hdbcNodes[i, 0] - 1, 0] - 0.5*apoio,
                 coords[hdbcNodes[i, 0] - 1, 1], marker='>', color='black',
                 markersize=8, markeredgecolor='blue')
    elif hdbcNodes[i, 1] == 2:
        ax1.plot(coords[hdbcNodes[i, 0] - 1, 0],
                 coords[hdbcNodes[i, 0] - 1, 1] - apoio, marker='^', color='black',
                 markersize=8, markeredgecolor='blue')

# Plot reaction forces
ax2 = fig1.add_subplot(212)
ax2.set_title('Reaction Loads', fontsize = 17)
plt.xlabel('X', fontsize = 10)
plt.ylabel('Y', fontsize = 10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax2)

# plot reaction loads
fFR = 0.25 * (rf * abs(max(prop[:, 0]))) / abs(max(rf))
for i in range(len(hdbcNodes)):
    ffr = np.asscalar(fFR[i])
    if hdbcNodes[i, 1] == 1:
        ax2.quiver(coords[hdbcNodes[i, 0] - 1, 0] - fFR[i, 0],
                   coords[hdbcNodes[i, 0] - 1, 1], fFR[i, 0], 0,
                   color='blue', scale=1, units='xy',
                   scale_units='xy', headlength=10, headwidth=5)
    elif hdbcNodes[i, 1] == 2:
        ax2.quiver(coords[hdbcNodes[i, 0] - 1, 0],
                   coords[hdbcNodes[i, 0] - 1, 1] - fFR[i, 0],
                   0, fFR[i, 0], color='blue', scale=1, units='xy',
                   scale_units='xy', headlength=10, headwidth=5)

nodes_mesh(ax2, order=totalNumElements + 1)
fig1.show()

### Plot Element Strains
fig2 = plt.figure()
plt.subplots_adjust(hspace=0.7, wspace=0.7)
ax3 = fig2.add_subplot(211)
ax3.set_title('Element Strains', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax3)

# Plot deformed mesh
# plot deformed mesh and element strains
d_ef = np.abs(coords - coordsd)
d_ef = np.sum(d_ef, axis=1)
z = rankmin(d_ef)
colors = pl.cm.jet(np.linspace(0, 1, len(set(z))))
Z = np.unique(z, return_index=False)
for elemNum in range(totalNumElements):
    xy = coordsdd[incid[elemNum, :] - 1, :]
    plot_gradient_rbg_pairs(xy[0, :], xy[1, :],
                            colors[np.where(z[incid[elemNum, 0] - 1] == Z)[0][0]],
                            colors[np.where(z[incid[elemNum, 1] - 1] == Z)[0][0]])

sm = plt.cm.ScalarMappable(cmap=plt.cm.jet,
                           norm=plt.Normalize(vmin=d_ef.min(), vmax=d_ef.max()))
sm._A = []
plt.colorbar(sm, ax=ax3, ticks=np.linspace(d_ef.min(), d_ef.max(), num=10))

# plot nodes of deformed mesh
nodes_mesh(ax3, coordsdd, order=totalNumElements + 1)

# Plot Element Stresses
ax4 = fig2.add_subplot(212)
ax4.set_title('Element stresses', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax4)

# plot stresses on elements of the deformed mesh
z2 = rankmin(stress)
colors2 = pl.cm.jet(np.linspace(0, 1, len(set(z2))))
Z2 = np.unique(z2, return_index=False)
for elemNum in range(totalNumElements):
        ax4.plot(coordsdd[incid[elemNum, :] - 1, 0],
                 coordsdd[incid[elemNum, :] - 1, 1],
                 linewidth=1, alpha=1, c=colors2[np.where(z2[elemNum] == Z2)[0][0]])

sm = plt.cm.ScalarMappable(cmap=plt.cm.jet,
                           norm=plt.Normalize(vmin=stress.min(), vmax=stress.max()))
sm._A = []
plt.colorbar(sm, ax=ax4, ticks=np.linspace(stress.min(), stress.max(), num=10))

# plot nodes of deformed mesh
nodes_mesh(ax4, coordsdd, order=totalNumElements + 1)
plt.show()

### Output File
outFileName = inpFileName + '.out2'
with open(outFileName, 'w') as f:
    f.write('*DISPLACEMENTS\n')
    for i in range(totalNumNodes):
        f.write('{:d} {:.4f} {:.4f}\n'.format(i + 1, u[nodalDOFNumbers[0, i] - 1, 0],
                                              u[nodalDOFNumbers[1, i] - 1, 0]))
    f.write('\n*ELEMENT_STRAINS\n')
    for i in range(totalNumElements):
        f.write('{:d} {:.6e}\n'.format(i + 1, strain[i, 0]))
    f.write('\n*ELEMENT_STRESSES\n')
    for i in range(totalNumElements):
        f.write('{:d} {:.6e}\n'.format(i + 1, stress[i, 0]))
    f.write('\n*REACTION_FORCES\n')
    for i in range(totalNumRestrDOFs):
        if hdbcNodes[i, 1] == 1:
            f.write('{:d} FX = {:.6e}\n'.format(hdbcNodes[i, 0], rf[i, 0]))
        else:
            f.write('{:d} FY = {:.6e}\n'.format(hdbcNodes[i, 0], rf[i, 0]))
    f.close()
