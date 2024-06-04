This is an example of two objective optimization using PVPmoo.
Here, the objective functions are final transverse projected emittance (mm-mrad)
and longitudinal bunch length (mm) from the ImpactT simulation of a 
photoinjector. The input file for the ImpactT simulation is copied
to the file "Input.in". There are five control parameters in this input
file will be adjuested to optimize the final two objective functions.
The locations of these five control parameters in the Input.in file are given
in the pvpmoo.in file:
------
152 157 163 163 171 /row line # of the control parameters in the input.in file
8 16 6 8 6 /column variable # of the control parameters in the input.in file
------
The range of these control parameters are given in the pvpmoo.in file:
------
20.0 0.02 -15.0d6 40.0 -9.0d6 /lower bounds of the control parameters
30.0 0.03 -10.0d6 50.0 -7.0d6 /upper bounds of the control parameters
------
Some other input files such as gun, solenoid, cavity field fields and the
ImpactT executable file (ImpactTexe) are also in this directory. 
During the optimization, a number of temporary directories (tmp0,...tmpNp-1)
will be created. Here, the number of directories is the same as the
number of processors used in the optimization. Here, the population in
the evolution algorithms is uniformly distributed among # of processors.
The ImpactT simulation is run in serial mode (i.e. on a single processor).
The optimization output is saved in obj.out and pop.out that contains nondominated
objective function values and control variable values of each generation.
The ImpactT input file ImpactT.in will be created based on the Input.in in
each temporary directory (tmp0,tmp1,...). The ImpactT simulation for each
injector configuration will be executed in the temporary directory. 

To run the optimizer, use "srun -n 16 ./xpvp" or "mpirun -n 16 ./xpvp".
Here, 16 is the number of processors used in the optimizer. 

It is OK to see the printout "cp: -r not specified; omitting directory 'tmp0'", etc.

