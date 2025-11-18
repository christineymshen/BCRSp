# slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# specnum <- as.integer(Sys.getenv("specnum"))

slurm_id <- 1
specnum <- 1

# folder path to store the simulation results
outfolder <- ""
# folder path which stores "functions.R" and the spec file
infolder <- ""

library(rstan)

runnum <- 1 # in case multiple runs were needed for the same spec
nsim <- 500 # number of simulations per slurm_id
isoverwrite <- T # whether overwrite existing saved run results
source(paste0(infolder,"functions.R"))

if (!file.exists(outfolder)) dir.create(outfolder)

spec <- readRDS(paste0(infolder, "spec_",specnum,".rds"))
n <- spec$n
start_seed <- (slurm_id-1)*nsim
mod <- stan_model(paste0(infolder, "CRS_is_HSGP5.stan"))

for (i in 1:nsim){
  seed <- start_seed+i
  output_name <- paste0(outfolder,"res",specnum,"_",runnum,"_",seed,".rds")
  
  if (!file.exists(output_name) | isoverwrite){
    set.seed(seed)
    sim <- simdata10(spec,seed=seed)
    max_time <- ceiling(max(sim$data$time))
    if (spec$pl){
      s <- sim$s[sim$s<=max_time]
      k <- length(s)-1
    } else {
      s <- seq(0,max_time,length.out=spec$k+1)
      k <- spec$k
    }
    
    runspec <- list(isModelrt=spec$isModelrt,isModeli=spec$isModeli,
                    isModels=spec$isModels,isJointwithin=spec$isJointwithin,
                    S=rep(spec$bdd,2),maxiter=30,y=sim$data$time,s=s,k=k,delta=sim$delta)
    
    j <- 1
    runspec <- updateHSGP_v2(runspec,spec,j)
    
    mod_time <- system.time({
      while (!runspec$check[[j]] & j<=runspec$maxiter){
        
        runspec <- runHSGP_v3(runspec,spec,mod,j,seed=seed)
        j <- j+1
        runspec <- updateHSGP_v2(runspec,spec,j)
        
      }
    })
    
    res <- runspec$fit[[j-1]]
    runspec$fit <- NULL
    output <- list(res=res,runspec=runspec,mod_time=mod_time)
    
    saveRDS(output, output_name)
  }
  
}
