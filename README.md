# PSOM_offline
This Fortran code advects tracers offline in PSOM model fields. 

***********
* INSTALLATION DIRECTIONS
This code is designed to work with pre-existing PSOM installations. In order to use this code, place the downloaded folder in your PSOM directory and compile using the same tools as your model run.

***********
* DATA NEEDS
Your namelist should include a path to the directory with the output of your model. These files at a minimum should include the face velocities (in 'face_xxxxx.cdf') and the free surface height (in 'full_xxxxx.cdf'). When doing offline advection, pay particular attention the output interval of the model (out3d_int). The code interpolates linearly between the input velocities.

***********
* NOTES ABOUT THE MODEL GRID
There are two ways that you set up the model grid for the offline advection. The first is to use the same files that you used to generate the model grid in the original model run. The second is to read the face coordinates into the model. In order to use this option, move the files in the subdirectory 'for_reading_grid' into the main directory.
