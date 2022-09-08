Script_SpecLocal produces simualtions and data analysis results.
Instructions for running analyses:
- Load R function files
- Alter paths in Script_SpecLocal to your local directory
- All write statements and figure pdf prints are commented out to avoid overwriting

There are some particular settings needed, which can be altered by commenting/uncommenting
sections of the script.  These are indicated in the comments of Script_SpecLocal.
- Adjust sample size T
- Determine which of four types of processes
- Determine parameters of the process

The script generates output for the basic and positive spectral estimator, based on the quadratic
as well as constant and quartic versions.  Various types of delta are allowed, as described in
the paper: an estimated delta, a theoretical delta (based on knowing the process), and a range of
delta values in increments of .005.  
