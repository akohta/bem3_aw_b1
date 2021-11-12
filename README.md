# bem3_aw_b1
This is the three-dimensional acoustic wave scattering analysis program for arbitrary objects in a inviscid fluid. 
This is based on boundary element method, the own developed numerical solution is used. 
The radiation force acting on the object can be analyzed. 
Intel Math Kernel Library and libpng are required. 
Gmsh is used for create a mesh data of object. 
The acoustic wave analysis program "multi_aw" is used for analyze incident field.

## Usage of example code  
1. type 'make' command to compile.  
   The executable aw_d3b1_bv_solver, example1.out, example2.out, example3.out are created. 
   The executable aw_d3b1_bv_solver is the main solver of boundary integral equations. 
   The example1.out is the executable of source code example1.c, it shows a simplest example using "bem3_aw_b1". 
   The example2.out is the execubable of source code example2.c, it shows a example of sound pressure intensity analysis. 
   The example3.out is the executable of source code example3.c, it shows a example of outputting the instantaneous value of sound pressure as an image.  
   
2. type './aw_d3b1_bv_solver' with arguments of medium datafile name, mesh datafile name and output dafafile name.  
   For example, './aw_d3b1_bv_solver medium.txt sphere_m1_05z.msh ex.dat'. 
   The medium.txt is the sample of medium datafile, one medium is defined in it. 
   The domain numbers are assigned to the medium from 1 in order. 
   The sphere_m1_05z.msh is the example of mesh datafile, it is a compressed sphere in the z direction. 
   It was created by Gmsh geometry file sphere_m1.geo in the mesh_sample folder and the stretch program of a object in the mesh_sample/stretch. 
   The sphere_m1_05z_image.png is the visulalization result of the sphere_m1.msh. 
   The aw_d3b1_bv_solver solves boundary integral equations with the specified datafiles, outputs the results to binary file with the output datafile name. 
   The mfb.txt is the sample of incident field datafile, a focused beam is defined in it. 
   Please refer to "multi_aw" for detail. 
   The aw_d3b1_bv_solver has optional arguments for rotation and translation of the object. 
   When the vector defining rotation axis is (rx, ry, rz), the rotation angle is theta, the translation vector is (tx, ty, tz), 
   the arguments are './aw_d3b1_bv_solver medium_data.txt sphere_m1.msh ex.dat rx ry rz theta tx ty tz' (these arguments are real number). 
   Rodrigues' rotation formula is used.  
   
3. type './example1.out' with an argument of datafile name output by aw_d3b1_bv_solver.  
   For example, './example1.out ex.dat'. 
   This executable calculates sound pressure, particle velocity, radiation force and torque.  
   
4. type './example2.out' with an argument of datafile name output by aw_d3b1_bv_solver.  
   For example, './example2.out ex.dat'. 
   This executable calculates sound pressure intensity distributions, outputs them to text files. 
   The I_example2.png is the visualization result of the intensity distributions, created by using Gnuplot script gscript_example2.plt.  
   
5. type './example3.out' with an argument of datafile name output by aw_d3b1_bv_solver.  
   For example, './example3.out ex.dat'. 
   This executable calculates instantaneous value of the sound pressure, outputs them to png image files. 
   The image files are output to the folder which has a name adding "images" to the datafile name specified in the argument (file-extension is excluded). 
   Each image file has a name that indicates the cross section, field component and number of time steps (ex. xz_p_014.png). 
   The color bar is output as color_bar.png in the same folder. 
   The range of color bar in each cross section is output to the info.txt file (ex. xy_info.txt for z=0 plane). 
   The xz_p.gif and xy_p.gif are animated gifs that concatenate the png files created by using the shell script gif_animation.sh.  
   
Please see d3b1_src/bem3_aw_b1.h for detail of functions. The main parts of the code are parallelized by using OpenMP. 
The number of threads is controlled by the environment variable OMP_NUM_THREADS. 
The additional analysis examples are in the folder analysis_sample1 ~ analysis_sample2.

![mesh 0](sphere_m1_05z_image.png "mesh image of the object (sphere_m1_05z_image.png)")  
![intensity distributions 0](I_example2.png "sound pressure intensity distributions (I_example2.png)")  
![xz_p.gif](xz_p.gif "instantaneous value of the p on y=0 plane (xz_p.gif)")![xy_p.gif](xy_p.gif "instantaneous value of the p on z=0 plane (xy_p.gif)")  


## Analysis sample 1 (in the folder analysis_sample1)  

This is the analysis result of focused beam scattering by the two layered sphere (silicone oil droplet containing air bubble).  

![mesh 1](analysis_sample1/sphere_m2_image.png "mesh image of the object (analysis_sample1/sphere_m2_image.png)")  
![intensity distributions 1](analysis_sample1/I_example2.png "sound pressure intensity distributions (analysis_sample1/I_example2.png)")  
![xz_p.gif 1](analysis_sample1/xz_p.gif "instantaneous value of the p on y=0 plane (analysis_sample1/xz_p.gif)")![xy_p.gif 1](analysis_sample1/xy_p.gif "instantaneous value of the p on z=0 plane (analysis_sample1/xy_p.gif)")  


## Verifications  
### Verification 1    
The first verification result using "aw_msp_ivf" is in the folder verification1.
It is the analysis result of focused beam scattering by the single sphere.
The sphere_m1_image.png is the visualization result of the mesh datafile. 
The I_example2.png is the visualization result of the intensity distributions.
The result of "aw_msp_ivf" is in the folder aw_msp_ivf_result.  

![mesh v1](verification1/sphere_m1_image.png "mesh image of the single sphere (verification1/sphere_m1_image.png)")  

### Verification 2  
The second verification result using "aw_msp_ivf" is in the folder verification2.
It is the analysis result of focused beam scattering by the three arranged spheres.
The sphere_m3_image.png is the visualization result of the mesh datafile. 
The I_example2.png is the visualization result of the intensity distributions.
The result of "aw_msp_ivf" is in the folder aw_msp_ivf_result.  

![mesh v2](verification2/sphere_m3_image.png "mesh image of the three arranged spheres (verification2/sphere_m3_image.png)")  


## About mesh file

This code can use quadrangular (bi-linear) and triangular (linear triangular) elements. 
I recommend using quadrangular element for reduce required memory. 
The samples of mesh datafile are in the folder mesh_sample. 
The file with extension .geo is the Gmsh geometry file. 
The file with extension .msh is the mesh datafile created by Gmsh geometry file. 
These mesh files are created by the command 'gmsh -2 -tol 1.0e-15 xxxx.geo' in command line (xxxx.geo is a geometry file). 
The domain number (Physical Surface) 99 is assigned to the open region in Gmsh geometry file, because Gmsh can't use the number 0 (assigned to open region in the code). 
Please refer to the manual of Gmsh for detail of geometry file.  


## References

1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)  
2. The official PNG reference library [libpng](http://www.libpng.org/pub/png/libpng.html)  
3. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)  
4. The command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
5. The utilities for manipulating images [ImageMagick](https://imagemagick.org/)  
6. The sound pressure analysis program [multi_aw](https://github.com/akohta/multi_aw/)  
7. The sound wave scattering analysis program [aw_msp_ivf](https://github.com/akohta/aw_msp_ivf/)  
