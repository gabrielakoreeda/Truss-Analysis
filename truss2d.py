### Importing libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.colors import LinearSegmentedColormap

### Reading the input filename
br = 60 * '='
print(br + '\n' + 'TRUSS 2D' + '\n' + br + '\n' + '\n' + br + '\n' + br)
InpFileName = 'example2_truss'
FileName = InpFileName + '.fem'
with open(FileName, 'r') as myFile:
    File = myFile.read().splitlines()
myFile.close()

### Reading function
def readNum(keyword, file, type = 'int'):
    if keyword in file:
        TotalNum = int(file[file.index(keyword) + 1])
        Num = pd.DataFrame(list(map(lambda x: x.split(), file[file.index(keyword) + 2: file.index(keyword) + 2 + TotalNum]))).astype(type).as_matrix()
        return TotalNum, Num
    else:
        raise ValueError('{} Not Found.'.format(keyword))
        
def removeFirstCol(dataframe):
    return dataframe[:, 1:]
    
### Reading nodal coordinates
print('READING NODAL COORDINATES...')
TotalNumNodes, Coords = readNum('*COORDINATES', File, type = 'float')
Coords = removeFirstCol(Coords)

# Number of nodal dofs
NumNodalDOFs = 2

### Reading element groups
print('READING ELEMENT GROUPS...')
NumGroups, Groups = readNum('*ELEMENT_GROUPS', File)
Groups = removeFirstCol(Groups)

# Total number of elements
TotalNumElements = Groups.sum()
    
### Reading incidences
print('READING ELEMENT INCIDENCES...')
if '*INCIDENCES' in File:
    Incid = pd.DataFrame(list(map(lambda groups: groups.split(), File[File.index('*INCIDENCES') + 1:File.index('*INCIDENCES') + TotalNumElements + 1]))).astype('int').as_matrix()
else:
    raise ValueError('{} Not Found'.format('*INCIDENCES'))
    
# deletes element numbers
Incid = removeFirstCol(Incid)
    
### Reading Materials
print('READING MATERIALS...')
NumMaters, Maters = readNum('*MATERIALS', File, type = 'float')

### Reading geometric properties
print('READING GEOMETRIC PROPERTIES...')
NumGPs, GPs = readNum('*GEOMETRIC_PROPERTIES', File, type = 'float')

### Reading boundary conditions
print('READING BOUNDARY CONDITIONS...')
NumBCNodes, HDBCNodes = readNum('*BCNODES', File)

### Reading loads
print('READING LOADS...')
NumLoadedNodes, Loads = readNum('*LOADS', File)
print(br)

### Element orientation: length, lx = sin theta, ly = cos theta
print('SOLUTION...\n' + br + '\n' + br)

# element lengths
lengths = np.sqrt(np.absolute((Coords[Incid[:, 0] - 1, 0] - Coords[Incid[:, 1] - 1, 0])**2 +
                   (Coords[Incid[:, 0] - 1, 1] - Coords[Incid[:, 1] - 1, 1])**2))
    
# lx = sin
lx = (Coords[Incid[:, 1] - 1, 1] - Coords[Incid[:, 0] - 1, 1])/lengths

# ly = cos
ly = (Coords[Incid[:, 1] - 1, 0] - Coords[Incid[:, 0] - 1, 0])/lengths

Prop = np.matrix([lengths, lx, ly]).transpose()

### Solver
# DOF Numbering
# DOFs matrix initialized with ones
NodalDOFNumbers = np.ones((NumNodalDOFs, TotalNumNodes), dtype = int)

# Boundary conditions on dofs with zero values 
NodalDOFNumbers[HDBCNodes[:NumBCNodes, 1] - 1, HDBCNodes[:NumBCNodes, 0] - 1] = 0

# Total number of DOFs for the model
TotalNumDOFs = TotalNumNodes*NumNodalDOFs

# Number of free and restricted dofs
TotalNumFreeDOFs = np.sum(NodalDOFNumbers == 1)
TotalNumRestrDOFs = np.sum(NodalDOFNumbers == 0)

# DOFs numbering
ukeqnum = 0
keqnum = TotalNumFreeDOFs
for i in range(TotalNumNodes):
    for j in range(NumNodalDOFs):
        if NodalDOFNumbers[j, i] == 1:
            ukeqnum = ukeqnum + 1
            NodalDOFNumbers[j, i] = ukeqnum
        elif NodalDOFNumbers[j, i] == 0:
            keqnum = keqnum + 1
            NodalDOFNumbers[j, i] = keqnum

### Assembling of the global stiffness matrix and load vector
# Creates zero global matrix
Kg = np.zeros((TotalNumDOFs, TotalNumDOFs))
ElemNum = 0
for Grp in range(NumGroups):
    E = Maters[Grp, 0]
    A = GPs[Grp, 0]
    for TNE in range(Groups[Grp - 1, 0]):
        cc = Prop[ElemNum, 2]**2
        ss = Prop[ElemNum, 1]**2
        cs = Prop[ElemNum, 1]*Prop[ElemNum, 2]
        
        # Element stiffness matrix
        Kedf = [[cc, cs, -cc, -cs],
                [cs, ss, -cs, -ss],
                [-cc, -cs, cc, cs],
                [-cs, -ss, cs, ss]]
        Ke = np.dot((E * A / Prop[ElemNum, 0]) , Kedf)
                
        # Element DOFs
        ElemEqs = NodalDOFNumbers[:, Incid[ElemNum, :] - 1]
        
        # Assembling into the global matrix
        index = np.concatenate(np.column_stack(ElemEqs)) - 1
        Kg[index, index.reshape((-1,1))] = Kg[index, index.reshape((-1,1))] + Ke
        ElemNum = ElemNum + 1
        
### Assembling global matrix
Fg = np.zeros((TotalNumDOFs, 1))  
for i in range(NumLoadedNodes):
    # element dofs
    ElemEqs = (NodalDOFNumbers[Loads[i, 1] - 1, Loads[i, 0] - 1])
    # assembling into the global vector
    Fg[ElemEqs - 1, 0] = Loads[i, 2]
    
### Solving systems of equations
U = np.zeros((TotalNumDOFs, 1))  
U[:TotalNumFreeDOFs, 0] = np.dot(np.linalg.inv(Kg[:TotalNumFreeDOFs, :TotalNumFreeDOFs]), Fg[:TotalNumFreeDOFs,0])
   
### Calculates reaction forces, element strains and stresses
# reaction forces
FR = np.dot(Kg[TotalNumFreeDOFs:, :TotalNumFreeDOFs], U[:TotalNumFreeDOFs])

# Deformed coordinates
Coordsd = Coords + U[np.concatenate(np.row_stack(NodalDOFNumbers)) - 1].reshape((NodalDOFNumbers.shape[0], NodalDOFNumbers.shape[1])).transpose()
ScaleFactor = 100
Coordsdd = Coords + ScaleFactor*U[np.concatenate(np.row_stack(NodalDOFNumbers)) - 1].reshape((NodalDOFNumbers.shape[0], NodalDOFNumbers.shape[1])).transpose()

# Element strain and stress vectors
Strain = np.zeros((TotalNumElements, 1))
Stress = np.zeros((TotalNumElements, 1))

# Calculates strain and stress for the elements
ElemNum = 0
for Grp in range(NumGroups):
    E = Maters[Grp, 0]
    for ElNum in range(Groups[Grp - 1, 0]):
        # element dofs and displacements
        ElemEqs = NodalDOFNumbers[:, Incid[ElemNum, :] - 1]
        index = np.concatenate(np.column_stack(ElemEqs)) - 1
        ElemU = U[index, 0]
        # strain
        Strain[ElemNum] = 1/Prop[ElemNum, 0]*np.dot(np.array([[-Prop[ElemNum, 2], -Prop[ElemNum, 1], Prop[ElemNum, 2], Prop[ElemNum,1]]]), ElemU)
        
        # stress
        Stress[ElemNum] = Strain[ElemNum]*E
        
        ElemNum = ElemNum + 1
  
### Post-Processing
print('POST-PROCESSING... \n' + br)

### Mesh and boundary conditions
plt.style.use('dark_background')
fig1 = plt.figure()
plt.subplots_adjust(hspace = 0.7, wspace = 0.7)

ax1 = fig1.add_subplot(211)
ax1.set_title('Truss', fontsize = 17)
plt.xlabel('X', fontsize = 10)
plt.ylabel('Y', fontsize = 10)
plt.xlim(min(Coords[:, 0]) - 0.5*Prop[0, 0], max(Coords[:, 0]) + 0.5*Prop[0, 0])
plt.ylim(min(Coords[:, 1]) - 0.5*Prop[0, 0], max(Coords[:, 1]) + 0.5*Prop[0, 0])
# plots original mesh
def originalMesh(ax):
    for ElemNum in range(TotalNumElements):
        ax.plot(Coords[Incid[ElemNum, :] - 1, 0], Coords[Incid[ElemNum, :] - 1, 1], color = 'white', linewidth = 1, alpha = 0.5)

originalMesh(ax1)
# plots nodes of the underformed mesh
def nodes_mesh(ax, coords = Coords):
    ax.scatter(x = coords[:, 0], y = coords[:, 1], c = 'yellow', s = 5, alpha = 0.5)

nodes_mesh(ax1)
### Plot loads
arrow_head = np.asscalar(0.4e-1*abs(max(Prop[:, 0])))
FFe = 0.25 * Loads[:, 2] * np.asscalar(abs(max(Prop[:, 0]))/abs(max(Loads[:, 2])))
for i in range(len(Loads)):
    if Loads[i, 1] == 1:
        ax1.arrow(Coords[Loads[i, 0] - 1, 0] - FFe[i], Coords[Loads[i, 0] - 1, 1], FFe[i], 0, color = 'red', head_width = arrow_head)
    elif Loads[i, 1] == 2:
        ax1.arrow(Coords[Loads[i, 0] - 1, 0], Coords[Loads[i, 0] - 1, 1] - FFe[i], 0, FFe[i], color = 'red', head_width = arrow_head)
        
### Plot bcs
apoio = np.asscalar(0.7e-1*abs(max(Prop[:, 0])))
for i in range(len(HDBCNodes)):
    if HDBCNodes[i, 1] == 1:
        ax1.arrow(Coords[HDBCNodes[i, 0] - 1, 0] - apoio, Coords[HDBCNodes[i, 0] - 1, 1], 1/100, 0, head_width = apoio, head_length = apoio, ec = 'blue', fc = 'none', lw = 2, )        
    elif HDBCNodes[i, 1] == 2:
        ax1.arrow(Coords[HDBCNodes[i, 0] - 1, 0], Coords[HDBCNodes[i, 0] - 1, 1] - apoio, 0, 0, head_width = apoio, head_length = apoio, ec = 'blue', fc = 'none', lw = 2)

### Plot reaction forces
ax2 = fig1.add_subplot(212)
ax2.set_title('Reaction Loads', fontsize = 17)
plt.xlabel('X', fontsize = 10)
plt.ylabel('Y', fontsize = 10)
plt.xlim(min(Coords[:, 0]) - 0.5*Prop[0, 0], max(Coords[:, 0]) + 0.5*Prop[0, 0])
plt.ylim(min(Coords[:, 1]) - 0.5*Prop[0, 0], max(Coords[:, 1]) + 0.5*Prop[0, 0])
originalMesh(ax2)  

# plot reaction loads
FFR = 0.25 * (FR * abs(max(Prop[:, 0]))) / abs(max(FR))
for i in range(len(HDBCNodes)):
    ffr = np.asscalar(FFR[i])
    if HDBCNodes[i, 1] == 1:
        ax2.arrow(Coords[HDBCNodes[i, 0] - 1, 0] - ffr, Coords[HDBCNodes[i, 0] - 1, 1], ffr, 0, head_width = arrow_head, color = 'blue')
    elif HDBCNodes[i, 1] == 2:
        ax2.arrow(Coords[HDBCNodes[i, 0] - 1, 0], Coords[HDBCNodes[i, 0] - 1, 1] - ffr, 0, ffr, head_width = arrow_head, color = 'blue')

nodes_mesh(ax2)   
fig1.show()     

# Color definition 
# n number of elements
def rankmin(x):
    x = np.round(x, decimals = 4)
    u, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    csum = np.zeros_like(counts)
    csum[1:] = counts[:-1].cumsum()
    return csum[inv]   

z = rankmin(Stress)
colors = pl.cm.jet(np.linspace(0, 1, len(set(z))))

### Plot Element Strains
fig2 = plt.figure()
plt.subplots_adjust(hspace = 0.7, wspace = 0.7)
ax3 = fig2.add_subplot(211)
ax3.set_title('Element Strains', fontsize = 17)  
plt.xlabel('X', fontsize = 10)
plt.ylabel('Y', fontsize = 10)
plt.xlim(min(Coords[:, 0]) - 0.5*Prop[0, 0], max(Coords[:, 0]) + 0.5*Prop[0, 0])
plt.ylim(min(Coords[:, 1]) - 0.5*Prop[0, 0], max(Coords[:, 1]) + 0.5*Prop[0, 0])
originalMesh(ax3)

def plot_gradient_hack( p0, p1, npts=20, cmap=None, **kw):
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

    C = np.linspace(0,1, npts)
    cmap = plt.get_cmap(cmap)
    # the simplest way of doing this is to just do the following:
    for x, y, c in zip(Xpairs, Ypairs, C):
        plt.plot(x, y, '-', c=cmap(c), **kw)

    # But for cases when that  will be too slow, you can make this go faster,
    # follow along with this example:
    # http://matplotlib.org/examples/pylab_examples/line_collection2.html


def plot_gradient_rbg_pairs(p0, p1, rgb0, rgb1, **kw):
    """Form the gradient from RGB values at each point
    The **kw dictionary gets passed to plt.plot, so things like linestyle,
    linewidth, labels, etc can be modified directly.
    """
    cmap = LinearSegmentedColormap.from_list('jet', (rgb0, rgb1))
    plot_gradient_hack(p0, p1, cmap=cmap, **kw)

# Plot deformed mesh
    # plot deformed mesh and element strains
d_ef = np.abs(Coords - Coordsd)
d_ef = np.sum(d_ef, axis = 1)
z2 = rankmin(d_ef)
colors2 = pl.cm.jet(np.linspace(0, 1, len(set(z2))))
Z2 = np.unique(z2, return_index = False)
for ElemNum in range(TotalNumElements):
    xy = Coordsdd[Incid[ElemNum, :] - 1, :]
    plot_gradient_rbg_pairs(xy[0, :], xy[1, :], colors2[np.where(z2[Incid[ElemNum, 0] - 1] == Z2)[0][0]], colors2[np.where(z2[Incid[ElemNum, 1] - 1] == Z2)[0][0]])

# plot nodes of deformed mesh
nodes_mesh(ax3, Coordsdd)  

### Plot Element Stresses
#z = Stress.ravel().argsort().argsort().reshape(Stress.shape)
ax4 = fig2.add_subplot(212)  
ax4.set_title('Element stresses', fontsize = 17)
plt.xlabel('X', fontsize = 10)
plt.ylabel('Y', fontsize = 10)
plt.xlim(min(Coords[:, 0]) - 0.5*Prop[0, 0], max(Coords[:, 0]) + 0.5*Prop[0, 0])
plt.ylim(min(Coords[:, 1]) - 0.5*Prop[0, 0], max(Coords[:, 1]) + 0.5*Prop[0, 0])
originalMesh(ax4)    
    
# plot stresses on elements of the deformed mesh
Z = np.unique(z, return_index = False)
for ElemNum in range(TotalNumElements):
        ax4.plot(Coordsdd[Incid[ElemNum, :] - 1, 0], Coordsdd[Incid[ElemNum, :] - 1, 1], linewidth = 1, alpha = 1, c = colors[np.where(z[ElemNum] == Z)[0][0]])

# plot nodes of deformed mesh
nodes_mesh(ax4, Coordsdd) 
plt.show()    

### Output File
OutFileName = InpFileName + '.out2'
with  open(OutFileName, 'w') as f:
    f.write('*DISPLACEMENTS\n')
    for i in range(TotalNumNodes):
        f.write('{:d} {:.4f} {:.4f}\n'.format(i + 1, U[NodalDOFNumbers[0, i] - 1, 0], U[NodalDOFNumbers[1, i] - 1, 0]))
    f.write('\n*ELEMENT_STRAINS\n')
    for i in range(TotalNumElements):
        f.write('{:d} {:.6e}\n'.format(i + 1, Strain[i, 0]))
    f.write('\n*ELEMENT_STRESSES\n')
    for i in range(TotalNumElements):
        f.write('{:d} {:.6e}\n'.format(i + 1, Stress[i, 0]))
    f.write('\n*REACTION_FORCES\n')
    for i in range(TotalNumRestrDOFs):
        if HDBCNodes[i, 1] == 1:
            f.write('{:d} FX = {:.6e}\n'.format(HDBCNodes[i, 0], FR[i, 0]))
        else:
            f.write('{:d} FY = {:.6e}\n'.format(HDBCNodes[i, 0], FR[i, 0]))
    f.close()
           
