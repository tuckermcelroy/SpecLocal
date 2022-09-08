# Format for output files

There are two types of simulations: (i) an extensive simulation involving four kinds of processes, (ii) a focused simulation of a Gaussian ARMA process with local quadratic, constant, and quartic estimators.

## Extensive simulation

There are 4 zip files containing output files for each type of process.
1. GaussianSims contains output files for Gaussian ARMA processes.
2. LaplaceSims contains output files for Laplace ARMA processes.
3. StudentSims contains output files for Student ARMA processes.
4. NonlinearSims contains output files for order 2 polynomial processes.

Each file has a label ending with T followed by either 50, 100, 200, 400, or 800, corresponding to the sample size.  There are two numbers before the T, 
each between 1 and 5, which indicate which specific process has been simulated.  For the first three kinds of processes, the numbers refer to the phi and theta
values in the ARMA specification.
- The first number is 1 corresponds to phi = -.9
- The first number is 2 corresponds to phi = -.5
- The first number is 3 corresponds to phi = 0
- The first number is 4 corresponds to phi = .5
- The first number is 5 corresponds to phi = .9
- The second number is 1 corresponds to theta = -.8
- The second number is 2 corresponds to theta = -.4
- The second number is 3 corresponds to theta = 0
- The second number is 4 corresponds to theta = .4
- The second number is 5 corresponds to theta = .8

Within each file there are 6 columns:
1. The first column corresponds to bias for frequency zero
2. The second column corresponds to standard deviation for frequency zero
3. The third column corresponds to square root mean squared error for frequency zero
4. The fourth column corresponds to bias for frequency pi
5. The fifth column corresponds to standard deviation for frequency pi
6. The sixth column corresponds to square root mean squared error for frequency pi

