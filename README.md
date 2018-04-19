# Truss-Analysis
This code written in Python analyses a 2D truss given an input file (.fem) that follows the following pattern.
*COORDINATES
6  
1 0.0 0.0
2 0.0 4.0
3 3.0 0.0
4 3.0 4.0
5 6.0 0.0
6 6.0 4.0

*ELEMENT_GROUPS
3
1 3
2 4
3 2
 
*INCIDENCES
1 1 2
2 3 4
3 5 6
4 1 3
5 3 5
6 2 4
7 4 6
8 1 4
9 5 4

*MATERIALS
3
210E9 120E6 80E6
80E9  70E6  60E6
70E9  85E6  56E6   

*GEOMETRIC_PROPERTIES
3
1E-4
2E-4
3E-4   

*BCNODES
5
1 1 
1 2
3 2
5 1 
5 2

*LOADS
4
2 2 -5000
4 2 -10000
6 2 -5000
6 1 3000

The code then calculates the strains, stress, reaction forces and displacements of each element of the truss and outputs these results in an output file like the following.
*DISPLACEMENTS
1 0.0000 0.0000
2 0.0010 -0.0010
3 0.0000 0.0000
4 0.0010 -0.0009
5 0.0000 0.0000
6 0.0014 -0.0013

*ELEMENT_STRAINS
1 -2.380952e-04
2 -2.352720e-04
3 -3.125000e-04
4 0.000000e+00
5 0.000000e+00
6 0.000000e+00
7 1.428571e-04
8 -3.152644e-05
9 -2.696217e-04

*ELEMENT_STRESSES
1 -5.000000e+07
2 -4.940711e+07
3 -2.500000e+07
4 0.000000e+00
5 0.000000e+00
6 0.000000e+00
7 1.000000e+07
8 -2.206851e+06
9 -1.887352e+07

*REACTION_FORCES
1 FX = 3.972332e+02
1 FY = 5.529644e+03
3 FY = 4.940711e+03
5 FX = -3.397233e+03
5 FY = 9.529644e+03

The Code also generates four plots:
1) the truss
2) the reaction forces
3) the strain
4) the stress
