# BCRSp

This repository contains reproducible codes for manuscript titled "Discovering Spatial Patterns of Readmission Risk Using a Bayesian Competing Risks Model with Spatially Varying Coefficients". In the manuscript, we described a simulation study and an application on Duke Electronic Health Record (EHR) data. Because the EHR data contain Protected Health Information (PHI), we cannot share these data. Here we have prepared:


1. In the `simulation` folder, `R` and `Stan` codes used for the simulation study. In our simulation study, we used the covariates and jittered locations of a subset of the EHR data to create the simulated datasets to ensure they closely resemble the real data. As we cannot share the EHR data, we have created a synthetic dataset based on the population characteristics of our EHR data, with locations uniformly drawn in a $2 \times 2$ box. We provide codes such that covariates and locations are extracted from this synthetic dataset to create competing risk datasets for simulation study. In this ways, the codes are executable. However, results will be different from those presented in our manuscript.
2. In the `application` folder, all the `R` and `Stan` codes we used for the application. These codes are not executable as we cannot make the underlying data available.

Descriptions on each of the files and how to use the codes to run simulation study are provided below. 

## Simulation

Here's a list of all the files and short descriptions:

1. `synthetic.rds`: Synthetic dataset which contains covariates and locations of $n=225$ observations. The covariates were generated based on the population characteristics of the subset of real EHR data we used for simulation study. The locations were uniformly chosen in the $[-1,1] \times [-1,1]$ box. With this file, the simulation codes can be run. However, as it's different from the real data that we used, the results will be different from what was presented in the manuscript.
2. `functions.R`: Customized functions used for the simulation study.
3. `CRS_is_GP2.stan`: `Stan` file for the proposed Bayesian competing risks spatial model using Gaussian Process (GP) priors.
4. `CRS_is_HSGP5.stan`: `Stan` file for the proposed Bayesisan competing risks spatial model using a Hilbert space low-rank approximation for GP (HSGP).
5. `1.sim_spec.R`: Codes to set up specifications for the simulation study. After executing this script, a `spec` file will be output to the user-specified folder.
6. `2.sim_freq.R`, `2.sim_GP.R` and `2.sim_HSGP.R`: Codes to generate 500 competing risks datasets using the synthetic data and then do model fitting, respectively using frequentist methods, full GP and HSGP.
7. `3.sim_summary.R`: Summarize model fitting results and produce figures presented in the manuscript.

Steps:

1. Save `synthetic.rds` and `functions.R` in the same folder, use this path as the input folder path in the few R script files, also use it as the output folder path in `1.sim_spec.R` file.
2. Run `1.sim_spec.R` file, a `spec1.rds` file will be output in the same folder.
3. Open each of the three `2.sim_[].R` files, specify an output folder for each one of them, and run. It's better to set up parallel runs for these simulation files, otherwise the two Bayesian runs will take a very long time.
4. After the run, update the relevant folder path in `3.sim_summary.R` and run this file to obtain summary figures as presented in the manuscript.

## Application

Here's a list of all the files and short descriptions:

1. `functions.R`: Customized functions used for the EHR data analysis.
2. `CRS_i_HSGP3.stan`: `Stan` file for the Bayesian competing risks model with spatial intercepts only.
3. `CRS_is_HSGP5.stan`: `Stan` file for the Bayesian competing risks model with both spatial intercepts and spatial slopes.
4. `CRS7.stan`: `Stan` file for the Bayesian competing risks model without spatial effects.
5. `1.spec_base.R": Create spec file for the base run.
6. `1.spec_highcorr.R`, `1.spec_lowcorr.R`, `1.spec_k100.R`: Create spec files for the sensitivity runs.
7. `2.run_nsp.R`: Script for Bayesian competing risks model run without spatial effects.
8. `2.run_sp.R`: Script for Bayesian competing risks model runs with spatial intercept, or with both spatial intercepts and slopes.
9. `3.summary_betatable.R`: Script to create Table 2 in the manuscript.
10. `3.summary_kmeans.R`: Script for k-mean clustering on the spatial random effects to create Figure 9 in the manuscript and related figures in the Supplement.
11. `3.summary_spatial.R`: Script to create Figure 8 in the manuscript and related figures in the Supplement.
12. `3.summary_waic_beta.R`: Script to produce cross validation results using WAIC, and 
