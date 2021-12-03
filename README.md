# KolmogorovFlowSwimmerDynamics
Code for calculating the steady-state distributions of noisy swimmers in a Kolmogorov flow, which accompanies the manuscript https://arxiv.org/abs/2111.09268.
Included is code for calculating the steady-state probability distribution of an ellipsoidal swimmer with rotational diffusion in a 2D laminar Kolmogorov flow, given the swimming speed, shape parameter, and rotational diffusivity of the swimmer.
In addition, there is code implementing the calculations of the averaging principle.
These functions average the noise-induced drift and diffusivity in phase space over the deterministic trajectories of the swimmer.
The averaged drift and diffusivities are then used to reconstruct the swimmer phase-space density in the weak-noise limit.

For a new user, download all files and directories. First, run plotKolmogorovPyth2.m, which generates data and makes the plot of Fig. 4 in the manuscript. Then, run plotKolmogorovAveragingCompare2.m, which uses data saved by the previous script, generates new data, and plot of Fig. 10 of the manuscript.
