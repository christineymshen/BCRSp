# BCRSp

This repository contains reproducible codes for manuscript titled "Discovering Spatial Patterns of Readmission Risk Using a Bayesian Competing Risks Model with Spatially Varying Coefficients". In the manuscript, we described a simulation study and an application on Duke Electronic Health Record (EHR) data. Because the EHR data contain Protected Health Information (PHI), we cannot share these data. Here we have prepared:


1. In the `simulation` folder, `R` and `Stan` codes used for the simulation study. In our simulation study, we used the covariates and jittered locations of a subset of the EHR data to create the simulated datasets to ensure they closely resemble the real data. As we cannot share the EHR data, we have created a synthetic dataset based on the population characteristics of our EHR data, with locations uniformly drawn in a $2 \times 2$ box. We provide codes such that covariates and locations are extracted from this synthetic dataset to create competing risk datasets for simulation study. In this ways, the codes are executable. However, results will be different from those presented in our manuscript.
2. In the `application` folder, all the `R` and `Stan` codes we used for the application. These codes are not executable as we cannot make the underlying data available.

Descriptions on each of the files and how to use the codes to run simulation study are provided below. 

## Simulation

Here's a list of all files and short descriptions:

1. `synthetic.rds`: Synthetic dataset which contains covariates and locations of $n=225$ observations. The covariates were generated based on the population characteristics of the subset of real EHR data we used for simulation study. The locations were uniformly chosen in the $[-1,1] \times [-1,1]$ box. With this file, the simulation codes can be run. However, as it's different from the real data that we used, the results will be different from what was presented in the manuscript.
2. `functions.R`: Contains customized functions used for the simulation study.
3. `CRS_is_GP2.stan`: `Stan` file used for the proposed Bayesian competing risks spatial model using Gaussian Process (GP) priors.
4. `CRS_is_HSGP5.stan`: `Stan` file used for the proposed Bayesisan competing risks spatial model using a Hilbert space low-rank approximation for GP (HSGP).
5. `1.sim_spec.R`: Codes to set up specifications for the simulation study. After executing this script, a `spec` file will be output to the user-specified folder.
6. `2.sim_freq.R`, `2.sim_GP.R` and `2.sim_HSGP.R`: Codes to generate 500 competing risks datasets using the synthetic data and then do model fitting, respectively using frequentist methods, full GP and HSGP.


## Application
