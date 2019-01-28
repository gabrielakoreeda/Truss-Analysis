# Truss-Analysis
This project applies the Finite Elements Method in the analysis of a truss and a portico type of structure.

## Introduction
There are three different codes written in Python. All of them require a input file of extension ".fem" which contains the properties of the structure. They also return a file ".out" that presents the  calculated strains, stress, reaction forces and displacements of each element.
Each code is separated in three parts: pre-processing, solver and post-processing.

### Pre-Processing
In this part, the necessary libraries are imported, the input file is read and the corresponding variables are attributed.

### Solver
The essential calculations are made to obtain the values of the strains, stress, reaction forces and displacements.

### Post-Processing
In the post-processing phase, plots are generated to better visualize the results. Also the output file is written.

## truss2d.py
This code analyses a simple truss under concentrated loads.

## truss2d_design.py
This code is a slight modification of the truss2d.py. It iterates over a maximum number of iterations defined to find the best geometric properties to bear the loads it is under.
## portico.py
This code analysis a truss that is under concentrated loads as well as distributed loads and momentum.
