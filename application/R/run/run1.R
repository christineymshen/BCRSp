specnum <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(rstan)
library(tidyverse)

runnum <- 1

folder <- "application/"

source(paste0(folder,"R/functions.R"))

spec <- readRDS(paste0(folder, "spec/spec",specnum,".rds"))

if (spec$isModels){
  stan_file <- paste0(folder,"stan/CRS_is_HSGP9.stan")
} else if (spec$isModeli) {
  stan_file <- paste0(folder,"stan/CRS_i_HSGP11.stan")
} else {
  stan_file <- paste0(folder,"stan/CRS7.stan")
}

stan_mod_file <- str_replace(stan_file, "\\.stan", "_mod.rds")
if (file.exists(stan_mod_file)){
  mod <- readRDS(stan_mod_file)
} else {
  mod <- stan_model(stan_file)
  saveRDS(mod, stan_mod_file)
}

if (spec$isModeli | spec$isModels){
  runspec <- list(S=spec$S,maxiter=30,y=spec$data$time,s=spec$s,delta=spec$delta,
                  nchain=4,niter=5000,nburn=1000,nthin=5)
  j <- 1
  runspec <- updateHSGP_v1(runspec,spec,j)
  
  mod_time <- system.time({
    while (!runspec$check[[j]] & j<=runspec$maxiter){
      
      runspec <- runHSGP_v5(runspec,spec,mod,j,seed=1)
      j <- j+1
      runspec <- updateHSGP_v1(runspec,spec,j)
      
    }
  })
  
  res <- runspec$fit[[j-1]]
  runspec$fit <- NULL
} else {
  runspec <- list(n=spec$n,p=spec$p,k=spec$k,m=spec$m,y=spec$data$time,
                  s=spec$s,delta=spec$delta,X=spec$X,
                  a0=spec$a0,b0=spec$b0,a1=spec$a1,b1=spec$b1,
                  nchain=4,niter=5000,nburn=1000,nthin=5)
  
  mod_time <- system.time(
  res <- sampling(mod, data=runspec, init=0.5, chains=runspec$nchain, 
                  warmup=runspec$nburn, iter=runspec$niter,core=runspec$nchain, 
                  thin=runspec$nthin, seed = specnum, refresh=0,
                  include=T, pars = c("lambda","log_lik","beta","kappa")) 
  )
}


output <- list(res=res,runspec=runspec,mod_time=mod_time)

saveRDS(output, paste0(folder,"res/spec",specnum,"_run",runnum,".rds")) 

