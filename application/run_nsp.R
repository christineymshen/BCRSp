start_specnum <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(rstan)
library(tidyverse)
nsim <- 4

runnum <- 1
folder <- ""
mod <- readRDS(paste0(folder,"stan/CRS7_mod.rds"))

source(paste0(folder,"R/functions.R"))

for (i in 1:nsim){
  specnum <- start_specnum + (i-1)*2
  spec <- readRDS(paste0(folder, "spec/spec",specnum,".rds"))
  
  if (spec$isModels){
    data <- list(n=spec$n,p=spec$p+1,k=spec$k,m=spec$m,y=spec$data$time,
                 s=spec$s,delta=spec$delta,X=cbind(spec$X,spec$W),
                 a0=spec$a0,b0=spec$b0,a1=spec$a1,b1=spec$b1)
  } else {
    data <- list(n=spec$n,p=spec$p,k=spec$k,m=spec$m,y=spec$data$time,
                 s=spec$s,delta=spec$delta,X=spec$X,
                 a0=spec$a0,b0=spec$b0,a1=spec$a1,b1=spec$b1)
  }
  
  res <- sampling(mod, data = data, init=0.5, chains = 4, warmup = 1000, 
                  iter = 5000, core=4, thin = 5, seed = specnum, 
                  include=T, pars = c("lambda","log_lik","beta","kappa")) 
  
  output <- list(res=res,runspec=data)
  
  saveRDS(output, paste0(folder,"res/CRS_np_spec",specnum,"_run",runnum,".rds"))
}
 
