specnum <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(rstan)
library(tidyverse)

runnum <- 1

wdir <- ""

if (!file.exists(wdir))
  wdir <- ""

folder <- paste0(wdir, "CRS/")

source(paste0(folder,"R/functions.R"))

spec <- readRDS(paste0(folder, "spec/spec",specnum,".rds"))

runspec <- list(S=spec$S,maxiter=30,y=spec$data$time,s=spec$s,delta=spec$delta,
                nchain=4,niter=5000,nburn=1000,nthin=5)

if (spec$isModels){
  mod <- readRDS(paste0(folder,"stan/CRS_is_HSGP5_mod.rds"))
} else {
  mod <- readRDS(paste0(folder,"stan/CRS_i_HSGP3_mod.rds"))
}

j <- 1
runspec <- updateHSGP_v1(runspec,spec,j)

mod_time <- system.time({
  while (!runspec$check[[j]] & j<=runspec$maxiter){
    
    runspec <- runHSGP_v2(runspec,spec,mod,j,seed=1,testfull=T)
    j <- j+1
    runspec <- updateHSGP_v1(runspec,spec,j)
    
  }
})

res <- runspec$fit[[j-1]]
runspec$fit <- NULL
output <- list(res=res,runspec=runspec,mod_time=mod_time)

saveRDS(output, paste0(folder,"res/HSGP_spec",specnum,"_run",runnum,".rds")) 

