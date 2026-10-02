# Bayesian competing risks model with spatially varying coefficients

Code accompanying the manuscript **“Discovering Spatial Patterns of Readmission Risk Using a Bayesian Competing Risks Model with Spatially Varying Coefficients.”**

The repository contains R and Stan code for the simulation study and the application to Duke electronic health record (EHR) data. The original EHR data contain protected health information and cannot be shared publicly. Therefore, we have prepared a synthetic dataset for illustration.

The synthetic dataset and two other data files used for visualization are provided in the `data` folder. Stan code is stored in the `stan` folder. In the `application` folder, we provide scripts used for the real-data application. Most of these scripts are for inspection only and cannot be executed using the public data alone. To illustrate the analysis, we also provide example scripts in the `application/R/example` folder. In the `simulation` folder, we provide scripts used for the simulation study. The simulation scripts use a subset of covariates and locations from the application data and generate competing-risks outcomes. Results obtained using synthetic data will differ from the manuscript results.

Detailed descriptions of the scripts and instructions for using them are provided below.

## Code structure

For both the application and the simulation, R scripts are organized into three subfolders: `run`, `spec`, and `summary`:

  - Model specification files (i.e., spec files) are in the `spec` folder. For example, spec 1 in the application folder is for our primary analysis (rather than the sensitivity analyses) with both spatial intercepts and spatial slopes, while spec 5 is for the primary analysis with spatial intercepts only. The mapping between spec numbers and models is recorded in the `speclog` file in the same folder. Running each of these scripts creates an `rds` file in the `application/spec` or `simulation/spec` folder. These `rds` files will later be used for model fitting and summarizing model results.

  - Model fitting files (i.e., run files) are in the `run` folder. These files will read in spec files based on the spec number specified by the user and then fit the model using Stan files. Run results will be stored in the `application/res` or `simulation/res` folders.

  - Results summary files (i.e., summary files) are in the `summary` folder. These files summarize and visualize run results for each model.

## Data

There are three files in the `data` folder.

  - `synthetic_v2.rds` contains a synthetic dataset with $n=1,200$ observations. The covariate distributions are similar to those in our real EHR data. Locations were uniformly drawn from the study region. This file is used in the application example scripts, as well as the simulation scripts.

  - The other two files contain geographical information that is used only to visualize model fitting results.

## Stan

There are four Stan files in this folder:

  - `CRS7.stan`: Bayesian competing risks model, without spatial random effects
  - `CRS_i_HSGP11.stan`: Bayesian competing risks model with spatially varying intercepts, implemented using an HSGP approximation to a GP
  - `CRS_is_HSGP9.stan`: Bayesian competing risks model with spatially varying intercepts and slopes, implemented using an HSGP approximation to a GP
  - `CRS_is_GP3.stan`: Bayesian competing risks model with spatially varying intercepts and slopes, implemented using a full GP

## Application

Scripts for the real-data application are stored in the `application/R` folder.

  - The `functions.R` file contains helper functions.
  - Scripts in the `run`, `spec`, and `summary` folders were used for the EHR data analysis. They are for inspection only and are not executable because we are not able to share the EHR data.
  - Scripts in the `example` folder are for illustrative purposes and can be executed. They are adapted from the scripts for our primary analysis using the proposed Bayesian model with both spatial intercepts and spatial slopes. The scripts were modified only to allow them to run on the synthetic data.

To use the example scripts:

1. Run the three `1.spec[specnum].R` files. `1.spec1.R` specifies the model for the primary analysis with both spatial intercepts and slopes, `1.spec5.R` specifies the model with only spatial intercepts, and `1.spec11.R` specifies the model without spatial random effects. After running these files, the corresponding model specification `rds` files will be saved in the `application/spec/` folder. Note that the `spec1.rds` file is already in the repository. This file can be created from the `1.spec1.R` file. The simulation specification scripts require this file, so we include it in case users try to run simulation scripts before using the application example scripts.

2. Run the `2.run1.R` file for model fitting. A total of three model runs are needed, one for each spec. Users can manually change the `specnum` parameter in this run file for different runs. The default setup uses four parallel chains, each with 5,000 iterations, including 1,000 warmup iterations. Users can adjust these parameters in the run file. Run results will be saved in the `application/res/` folder.

3. Run the `3.spec1_base.R` file to analyse and visualize the run results. Figures will be saved in the `application/fig/` folder. Note that because the synthetic data differ from the real EHR data, the results will be different from those presented in the manuscript.

## Simulation

Scripts for the simulation study are stored in the `simulation/R` folder.

  - The `functions.R` file contains helper functions.
  - Scripts in the `run`, `spec`, and `summary` folders were adapted from those used for our simulation study. The scripts were modified only to allow them to run directly on the synthetic data.

To run the simulation scripts:

1. Run the model specification files in the `spec` folder, and the corresponding `spec[specnum].rds` files will be created and saved in the `simulation/spec` folder. `spec1.R` is for a simulation study with sample size $n=225$. Most of the simulation results presented in the manuscript used this setup. `spec2.R` is for a sample size of $n=500$, and `spec3.R` is for a sample size of $n=100$. Simulation results under these two sample sizes were presented in the Supplementary Materials, Tables 2-4. Note that these spec files use the synthetic data indirectly through the `application/spec/spec1.rds` file.

2. Model fitting can be done using files in the `run` folder. There are three files in this folder: one (`run1_HSGP.R`) is for the Bayesian model using an HSGP approximation, one (`run1_GP.R`) is for the Bayesian model using a full GP, and one (`run1_FQ.R`) is for the two frequentist runs. For model fitting on the three different model specifications, users need to manually update the `specnum` parameter in these run files. We used SLURM array jobs on a computing server to implement these simulations. The run results will be saved in the `simulation/res` folder.

3. Note that the frequentist run needs to precede the Bayesian model runs. This is because when the sample size is small, there might be no events in the minority groups in some of the simulated datasets. The frequentist run might fail to converge for these datasets. To keep the results comparable, we will not use datasets that result in fitting errors, and will instead generate another dataset using a different seed. Only after finishing the frequentist runs do we know which seeds can be used for the simulated datasets, and the same seeds will be used for the Bayesian runs for comparability. During our implementation, this was only a problem for the $n=100$ set, and there were fewer than five seeds that we needed to discard due to this issue.

4. Use the R script files in the `summary` folder to summarize and visualize the simulation results.
