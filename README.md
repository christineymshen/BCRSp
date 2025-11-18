# BCRSp

This repository contains reproducible codes for manuscript titled "Discovering Spatial Patterns of Readmission Risk Using a Bayesian Competing Risks Model with Spatially Varying Coefficients". In the manuscript, we described a simulation study and an application on Duke Electronic Health Record (EHR) data. Because the EHR data contain Protected Health Information (PHI), we cannot share these data. Here we have prepared:


1. In the `simulation` folder, `R` and `Stan` codes used for the simulation study. In our simulation study, we used the covariates and jittered locations of a subset of the EHR data to create the simulated datasets to ensure they closely resemble the real data. As we cannot share the EHR data, we have created a synthetic dataset based on the population characteristics of our EHR data, with locations uniformly drawn in a $2 \times 2$ box. We provide codes such that covariates and locations are extracted from this synthetic dataset to create competing risk datasets for simulation study. In this ways, the codes are executable. However, results will be different from those presented in our manuscript.
2. In the `application` folder, all the `R` and `Stan` codes we used for the application. These codes are not executable as we cannot make the underlying data available.

Descriptions on each of the files and how to use the codes to run simulation study are provided below. 

## Simulation




## Application
