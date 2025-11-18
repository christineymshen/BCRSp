# BCRSp

This repository contains reproducible codes for manuscript titled "Discovering Spatial Patterns of Readmission Risk Using a Bayesian Competing Risks Model with Spatially Varying Coefficients". In the manuscript, we described a simulation study and an application on Duke Electronic Health Record (EHR) data. Because the EHR data contain Protected Health Information (PHI), we cannot share these data. Therefore here we provide:


1. In the `simulation` folder, `R` and `Stan` codes used for the simulation study. In our simulation study, we used the covariates and jittered locations of a subset of the EHR data to create the simulated datasets to ensure these datasets closely resemble the real data. As we cannot share the EHR data, we have created a synthetic dataset based on the population characteristics of our EHR data, with locations uniformly drawn in a $2 \times 2$ box. In this ways, the codes are executable. However, results will be different from those presented in our manuscript.
2. In the `application` folder, all the `R` and `Stan` codes we used for the application. These codes are not executable 
