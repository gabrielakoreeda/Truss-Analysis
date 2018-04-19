# Truss-Analysis
This code written in Python analyses a 2D truss given an input file (.fem) that follows the following pattern and just like the example upload.
*COORDINATES  

*ELEMENT_GROUPS
 
*INCIDENCES

*MATERIALS

*GEOMETRIC_PROPERTIES

*BCNODES

*LOADS

The code then calculates the strains, stress, reaction forces and displacements of each element of the truss and outputs these results in an output file like the following.
*DISPLACEMENTS

*ELEMENT_STRAINS

*ELEMENT_STRESSES

*REACTION_FORCES

The Code also generates four plots:
1) the truss
2) the reaction forces
3) the strain
4) the stress
