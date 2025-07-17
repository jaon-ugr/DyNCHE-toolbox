# DyNCHE-toolbox
DyNCHE-toolbox is dedicated to numerical resolution of the dynamics of a field on a box with moving boundaries. 

exec.c file is the main file where we run the program. It calls files in Library folder, which contain functions, odes, and initial data.

In functions.c we define the trajectories of the boundaries. The are written in a way to maximize efficiency. Hence, in order to change from one trajectory to another we chose to comment and uncomment lines. This must also be done in the odes.c file, since there one selects the concrete trajectory function for the boundaries.  

In order to run the code, we need gcc compiler and the gsl library. 

We then must give the script "script.sh" executable privileges. For instance

~ chmod a+x script.sh

It contains the names of the output files, the interval of modes that we will simulate, and the concrete initial mode "in solution" to be evolved. The last two parameters are the initial and final times. However, these parameters are currently redefined in the exec.c file, lines 42 and 43 for the oscillatory trajectories of the Dynamical Casimir effect and lines 46 and 47 for the Fulling and Davies trajectory (expanding cavity). 

Then, we run the script as

~ ./script.sh

We get output files "alpha-n16-1.txt" "beta-n16-1.txt" and "norm-n16-1.txt" where it prints the alpha and beta Bogoliubov coefficients and the norm of the solutions at the end for 16 total modes and initial solution peaked at "in" mode number 1. For other "in" modes we can change the last number in the "script.sh" file, choosing any number between 1 to 16. We can change the total number of modes by changing "16" to some other number.

