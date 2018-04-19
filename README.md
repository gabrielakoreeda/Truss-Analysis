# Truss-Analysis
This code written in Python analyses a 2D truss given an input file (.fem) that follows the following pattern and just like the example upload.

*COORDINATES  
6  

1 0.0 0.0

2 0.0 4.0

3 3.0 0.0

4 3.0 4.0

5 6.0 0.0

6 6.0 4.0

*ELEMENT_GROUPS
 
*INCIDENCES

*MATERIALS

*GEOMETRIC_PROPERTIES

*BCNODES

*LOADS

The code then calculates the strains, stress, reaction forces and displacements of each element of the truss and outputs these results in an output file like the following.

*DISPLACEMENTS

1 0.0000 0.0000

2 0.0010 -0.0010

3 0.0000 0.0000

4 0.0010 -0.0009

5 0.0000 0.0000

6 0.0014 -0.0013

*ELEMENT_STRAINS

*ELEMENT_STRESSES

*REACTION_FORCES

The Code also generates four plots:
1) the truss
2) the reaction forces
3) the strain
4) the stress
