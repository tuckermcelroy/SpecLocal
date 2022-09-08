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

For the nonlinear process, the numbers refer to the two values of phi in the specification.
- The first number is 1 corresponds to first phi = -.9
- The first number is 2 corresponds to first phi = -.5
- The first number is 3 corresponds to first phi = 0
- The first number is 4 corresponds to first phi = .5
- The first number is 5 corresponds to first phi = .9
- The second number is 1 corresponds to second phi = -.9
- The second number is 2 corresponds to second phi = -.5
- The second number is 3 corresponds to second phi = 0
- The second number is 4 corresponds to second phi = .5
- The second number is 5 corresponds to second phi = .9

Within each file there are 6 columns:
1. The first column corresponds to bias for frequency zero
2. The second column corresponds to standard deviation for frequency zero
3. The third column corresponds to square root mean squared error for frequency zero
4. The fourth column corresponds to bias for frequency pi
5. The fifth column corresponds to standard deviation for frequency pi
6. The sixth column corresponds to square root mean squared error for frequency pi

Within each file there are 106 rows:
- Row 1 corresponds to the Parzen taper
- Row 2 corresponds to the flat-top taper
- Row 3 corresponds to local quadratic with estimated delta
- Row 4 corresponds to local quadratic with theoretical delta
- Rows 5-54 correspond to local quadratic with delta in range .005, .010, ..., .250
- Row 55 corresponds to positive local quadratic with estimated delta
- Row 56 corresponds to positive local quadratic with theoretical delta
- Rows 57-106 correspond to positive local quadratic with delta in range .005, .010, ..., .250

## Focused simulation

The file labels correspond to three sample sizes (50, 200, 800) and the 54 ARMA Gaussian process, which has phi = .9 and theta = .8.

Within each file there are 6 columns:
1. The first column corresponds to bias for frequency zero
2. The second column corresponds to standard deviation for frequency zero
3. The third column corresponds to square root mean squared error for frequency zero
4. The fourth column corresponds to bias for frequency pi
5. The fifth column corresponds to standard deviation for frequency pi
6. The sixth column corresponds to square root mean squared error for frequency pi

Within each file there are 314 rows:
- Row 1 corresponds to the Parzen taper
- Row 2 corresponds to the flat-top taper
- Row 3 corresponds to local constant with estimated delta
- Row 4 corresponds to local constant with theoretical delta
- Rows 5-54 correspond to local constant with delta in range .005, .010, ..., .250
- Row 55 corresponds to positive local constant with estimated delta
- Row 56 corresponds to positive local constant with theoretical delta
- Rows 57-106 correspond to positive local constant with delta in range .005, .010, ..., .250
- Row 107 corresponds to local quadratic with estimated delta
- Row 108 corresponds to local quadratic with theoretical delta
- Rows 109-158 correspond to local quadratic with delta in range .005, .010, ..., .250
- Row 159 corresponds to positive local quadratic with estimated delta
- Row 160 corresponds to positive local quadratic with theoretical delta
- Rows 161-210 correspond to positive local quadratic with delta in range .005, .010, ..., .250
- Row 211 corresponds to local quartic with estimated delta
- Row 212 corresponds to local quartic with theoretical delta
- Rows 213-262 correspond to local quartic with delta in range .005, .010, ..., .250
- Row 263 corresponds to positive local quartic with estimated delta
- Row 264 corresponds to positive local quartic with theoretical delta
- Rows 265-314 correspond to positive local quartic with delta in range .005, .010, ..., .250
