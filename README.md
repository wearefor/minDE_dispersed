# COMSOL Simulation of Dispersed-Membrane MinDE Model
This repository contains the files necessary to recreate the simulations in "3D pattern formation of a protein-membrane suspension" (https://arxiv.org/abs/2501.03179). 

Experimental work by Amelie Chardac (a) and Guillaume Duclos (a)

Theory and computational work by Michael M. Norton (a) and Jonathan Touboul (b)

(a) Martin A. Fisher School of Physics, Brandeis University, Waltham, MA, USA

(b) Department of Mathematics, Brandeis University, Waltham, MA, USA



# I. COMSOL *.mph

The file "emptycomsol.mph" is an empty COMSOL simulation file, ready to run. Adjust parameters {nD,nE,alpha,beta,gamma,etc..} as desired. In the event that this file is incompatible with your version of COMSOL, the MATLAB file below can possibly be used to generate an initial *.mph file, since the core commands are highly compatible across COMSOL versions.

# II. MATLAB Livelink

"matlabsweep.m" will run a batch of simulations in serial. It is configured to run a sweep over a two-dimensional grid of nD and nE values.
For each simulation, a set of images and an animated *.gif of the total minD concentration will be produced along with two *.mph files (an empty one and one with data). 
Steps:
1. Initiate the COMSOL server using the following command
2. From within MATLAB, run the command "mphstart" or, if the port number must be specified, "mphstart(2036)". If this step fails, ensure that COMSOL's livelink library has been added to MATLAB's path. 


