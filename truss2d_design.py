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
                              )
                         ).astype('int').as_matrix()
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

# Reading maximum design iterations
print('READING DESIGN ITERATIONS...')
if '*DESIGN_ITERATIONS' in file:
    maxNumDesignIter = int(file[file.index('*DESIGN_ITERATIONS') + 1])
else:
    raise ValueError('{} Not Found.'.format('*DESIGN_ITERATIONS'))

# Element orientation: length, lx = sin theta, ly = cos theta
print('SOLUTION...\n' + br + '\n' + br)

# element lengths
lengths = np.sqrt(np.absolute((coords[incid[:, 0] - 1, 0] -
                               coords[incid[:, 1] - 1, 0])**2 +
                              (coords[incid[:, 0] - 1, 1] -
                               coords[incid[:, 1] - 1, 1])**2))
    
# lx = sin
lx = (coords[incid[:, 1] - 1, 1] - coords[incid[:, 0] - 1, 1]) / lengths

# ly = cos
ly = (coords[incid[:, 1] - 1, 0] - coords[incid[:, 0] - 1, 0]) / lengths

prop = np.matrix([lengths, lx, ly]).transpose()

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

# Assembling of the global stiffness matrix and load vector
# Areas and Materials Properties for each element
areas = np.zeros((totalNumElements, maxNumDesignIter))
matProp = np.zeros((totalNumElements, 3))
elemNum = 0
for grp in range(numGroups):
    for el in range(int(groups[grp])):
        areas[elemNum, 0] = gps[grp, 0]
        matProp[elemNum, :] = maters[grp, :3]
        elemNum = elemNum + 1

# Element Strain and Stress Vectors
strain = np.zeros((totalNumElements, maxNumDesignIter))
stress = np.zeros((totalNumElements, maxNumDesignIter))

# Assembling global load vector
fg = np.zeros((totalNumDOFs, 1))
for i in range(numLoadedNodes):
    # element dofs
    elemEqs = nodalDOFNumbers[loads[i, 1] - 1, loads[i, 0] - 1]
    
    # assembling into the global vector
    fg[elemEqs - 1, 0] = loads[i, 2]
    
# Global displacement vector
u = np.zeros((totalNumDOFs, 1))

# Sets number of design iterations and design flag (=1, new design; 0= final design)
nDIt = 1
designFlag = 1

while designFlag:
    # Creates zero global matrix and load vector
    kg = np.zeros((totalNumDOFs, totalNumDOFs))
    for elemNum in range(totalNumElements):
        # Youngs modulus and cross section area
        e = matProp[elemNum, 0]
        a = areas[elemNum, nDIt - 1]
        
        # director co-sine
        cc = prop[elemNum, 2]**2
        ss = prop[elemNum, 1]**2
        cs = prop[elemNum, 1]*prop[elemNum, 2]
        
        # element stiffness matrix
        kedf = [[cc, cs, -cc, -cs],
                [cs, ss, -cs, -ss],
                [-cc, -cs, cc, cs],
                [-cs, -ss, cs, ss]]
        ke = np.dot((e * a / prop[elemNum, 0]), kedf)
        
        # element dofs
        elemEqs = nodalDOFNumbers[:, incid[elemNum, :] - 1]
        
        # assembling into the global matrix
        index = np.concatenate(np.column_stack(elemEqs)) - 1
        kg[index, index.reshape((-1, 1))] = kg[index, index.reshape((-1, 1))] + ke
        
    # Solving systems of equations
    u[:totalNumFreeDOFs, 0] = np.dot(np.linalg.inv(kg[:totalNumFreeDOFs,
                                                   :totalNumFreeDOFs]),
                                     fg[:totalNumFreeDOFs, 0])
    
    # Calculates strain and stress or the elements
    for elemNum in range(totalNumElements):
        
        # Youngs modulus
        e = matProp[elemNum, 0]
        
        # Elements DOFs and displacements
        elemEqs = nodalDOFNumbers[:, incid[elemNum, :] - 1]
        index = np.concatenate(np.column_stack(elemEqs)) - 1
        elemU = u[index, 0]
        
        # Strain
        strain[elemNum, nDIt - 1] = 1 / prop[elemNum, 0]*np.dot(np.array([[-prop[elemNum, 2],
                                                                           -prop[elemNum, 1],
                                                                           prop[elemNum, 2],
                                                                           prop[elemNum, 1]]]),
                                                                elemU)
        
        # Stress
        stress[elemNum, nDIt - 1] = strain[elemNum, nDIt - 1]*e
    
    # if the maximum number of design iterations is reached, finish the design process
    # checks if the stresses for all bars are below the admissable limits or
    # if the maximum number of iterations has been reached
    tractionBars = stress[:, nDIt - 1] >= 0
    compressionBars = stress[:, nDIt - 1] < 0
    if (sum(stress[tractionBars, 0] <= matProp[tractionBars, 1]) +
        sum(stress[compressionBars, 0] >= -matProp[compressionBars, 2]) == totalNumElements) \
            or (nDIt + 1 >= maxNumDesignIter):
            designFlag = 0
    # otherwise calculate the new areas
    else:
        # increments number of design iterations
        nDIt += 1
       
        # copies areas from the previous iteration
        areas[:, nDIt - 1] = areas[:, nDIt - 2]
       
        # bars under compression with stress less than admissable
        # compression stress. They will be redesign.
        newBars = stress[:, nDIt - 2] < - matProp[:, 2]
       
        # new areas for bars under compression
        areas[newBars, nDIt - 1] = (- stress[newBars, nDIt - 2] /
                                    matProp[newBars, 2] *
                                    areas[newBars, nDIt - 2])
       
        # bars under traction with stress larger than admissable
        # tensile stress. They will be redesigned.
        newBars = stress[:, nDIt - 2] > matProp[:, 1]
       
        # new areas for bars under traction
        areas[newBars, nDIt - 1] = (stress[newBars, nDIt - 2] /
                                    matProp[newBars, 1] *
                                    areas[newBars, nDIt - 2])
          
# reaction forces
rf = np.dot(kg[totalNumFreeDOFs:, :totalNumFreeDOFs], u[:totalNumFreeDOFs])

# Deformed coordinates
coordsd = coords + u[np.concatenate(np.row_stack(nodalDOFNumbers)) - 1].reshape(
    (nodalDOFNumbers.shape[0], nodalDOFNumbers.shape[1])).transpose()
scaleFactor = 100
coordsdd = coords + scaleFactor*u[np.concatenate(np.row_stack(nodalDOFNumbers)) - 1].reshape(
    (nodalDOFNumbers.shape[0], nodalDOFNumbers.shape[1])).transpose()

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

# Plots Nodes of the Underformed Mesh
nodes_mesh(ax1)

# Plot Loads
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
        
# Plot BCs
apoio = np.asscalar(0.7e-1*abs(max(prop[:, 0])))
for i in range(len(hdbcNodes)):
    if hdbcNodes[i, 1] == 1:
        ax1.plot(coords[hdbcNodes[i, 0] - 1, 0] - apoio,
                 coords[hdbcNodes[i, 0] - 1, 1], marker='>', color='black',
                 markersize=8, markeredgecolor='blue')
    elif hdbcNodes[i, 1] == 2:
        ax1.plot(coords[hdbcNodes[i, 0] - 1, 0],
                 coords[hdbcNodes[i, 0] - 1, 1] - apoio,
                 marker='^', color='black', markersize=8, markeredgecolor='blue')

# Plot reaction forces
ax2 = fig1.add_subplot(212)
ax2.set_title('Reaction Loads', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax2)  

# plot reaction loads
fFR = 0.25 * (rf * abs(max(prop[:, 0]))) / abs(max(rf))
for i in range(len(hdbcNodes)):
    if hdbcNodes[i, 1] == 1:
        ax2.quiver(coords[hdbcNodes[i, 0] - 1, 0] - fFR[i, 0],
                   coords[hdbcNodes[i, 0] - 1, 1], fFR[i, 0], 0,
                   color='blue', scale=1, units='xy',
                   scale_units='xy', headlength=10, headwidth=5)
    elif hdbcNodes[i, 1] == 2:
        ax2.quiver(coords[hdbcNodes[i, 0] - 1, 0],
                   coords[hdbcNodes[i, 0] - 1, 1] - fFR[i, 0], 0,
                   fFR[i, 0], color='blue', scale=1, units='xy',
                   scale_units='xy', headlength=10, headwidth=5)

nodes_mesh(ax2, order=totalNumElements + 1)
fig1.show()     

# Plot Element Strains
fig2 = plt.figure()
plt.subplots_adjust(hspace=0.7, wspace=0.7)
ax3 = fig2.add_subplot(211)
ax3.set_title('Element Strains', fontsize=17)
plt.xlabel('X', fontsize=10)
plt.ylabel('Y', fontsize=10)
plt.xlim(min(coords[:, 0]) - 0.5*prop[0, 0], max(coords[:, 0]) + 0.5*prop[0, 0])
plt.ylim(min(coords[:, 1]) - 0.5*prop[0, 0], max(coords[:, 1]) + 0.5*prop[0, 0])
originalMesh(ax3)

# Plot Deformed Mesh
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

# Output File
outFileName = inpFileName + 'py.out2'
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
    f.write('\n*AREAS\n')
    f.write('{:d}\n'.format(nDIt))
    for i in range(totalNumElements):
        f.write('{:d} '.format(i + 1))
        for j in range(nDIt):
            f.write('{:e} '.format(areas[i, j]))
        f.write('\n')
    f.write('\n*VOLUMES\n')
    f.write('{:d}\n'.format(nDIt))
    for i in range(nDIt):
        f.write('{:e}\n'.format(np.asscalar(np.dot(areas[:, i], prop[:, 0]))))
    f.close()
    
if nDIt >= maxNumDesignIter:
    print('\nThe maximum number of design iterations has been reached.\n '
          'The stress for some bars may be above the given admissible stresses!!\n')
