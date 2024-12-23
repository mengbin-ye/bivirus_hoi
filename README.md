# bivirus_hoi

This Git repository provides some simple code in MATLAB to simulate the Bivirus SIS network model with Higher-Order Interactions (HOIs). A relevant paper (pre-print) with associated notation is https://hal.science/hal-04669193/

There are two files.

bivirus_dynamics_sicon.m is the main run file, and running this will allow the bivirus model to be simulated for a range of different parameters. 

In the front part of the run file, one can set the various parameters associated with the bivirus model, including the number of nodes, the infection/recovery rates, as well as the infection matrices for the pairwise and higher-order interactions. The parameter selection is done by using the "case" environment, so that you need to set the switch number in order to select the particular caes. This allows retaining of different simulation scenarios, as well as rapid switching between scenarios for comparison or testing.

The second section uses ODE45 to simulate the bivirus dynamics up to a specific simulation time window that can be adjusted. The final equilibrium values are recorded, as are some potentially relevant quantities.

The final section provides a simple plot of the trajectories over time.


bivirus_hoi.m is the function file that contains the ODE for the bivirus network model. Changes here will change the actual dynamical system that is being simulated
