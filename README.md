# isosurface plotting
The package is used to plot snowflake crystal particles (from OpenSSP) by using the MATLAB function isosurface.

# File list
run_att.sh     :  run script for Linux/MacOS

run_att.bat    :  run script for Windows

process_att.m  :  main MATLAB code for plotting particle voxel file (.att from OpenSSP)

readAtt.m      :  reading .att file

stl_write.m    :  writing .stl file based on the read-in .att file

list_att_files.txt:  input file for process_att.m, with a list of .att file for plotting

p-08: a folder contains several sample particle files from OpenSSP

plots: a folder contanis several generated plots

# Introduction
A snowflake crystal particle file (./p-08/.att) contains a regular 3-D grid array in Cartesian Geometry, with each element representing the material type/density of the corresponding voxel (e.g. 0 for air, 1 for ice).


After reading in the 3-D array from a .att file, process_att.m performs four operations before rendering the plot:

1. morphological operation : (MATLAB imclose)

2. convolution with 3x3x3 cube and Gaussian weights based on distance to the cube center

3. generation of isosurface (MATLAB isosurface)  

4. Coplannar operation.            

These operations are intened to smooth the particle surface while retaining key details of the particle shape. The coplannar operation identifies groups of neighboring coplannar triangles by evaluating each triangle's normal vector. If the inner product of two normal verctors is greater than a user-specifed threshold (e.g. 0.999), the two corresponding triangles are considered to be in the same coplannar group. Only the outer boundary of a coplannar group of triangles will be plotted. 
   

