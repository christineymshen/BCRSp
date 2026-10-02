slurm_id <- "to be updated"
specnum <- "to be updated"
runnum <- 1

FQrunnum <- 1

isoverwrite <- F
runlabel <- "HSGP"

folder <- "simulation/"
source(paste0(folder,"R/functions.R"))
outfolder <- paste0(folder,"res/", runlabel)
if (!file.exists(outfolder)) dir.create(outfolder, recursive = T)

spec <- readRDS(paste0(folder, "spec/spec",specnum,".rds"))

stan_file <- "stan/CRS_is_HSGP9.stan"
stan_mod_file <- str_replace(stan_file, "\\.stan", "_mod.rds")
if (file.exists(stan_mod_file)){
  mod <- readRDS(stan_mod_file)
} else {
  mod <- stan_model(stan_file)
  saveRDS(mod, stan_mod_file)
}

output_name <- paste0(outfolder,"/spec",specnum,"_run",runnum,"_",slurm_id,".rds")

seed_idx <- readRDS(paste0(folder,"res/FQ/spec",specnum,"_run",FQrunnum,"_seed.rds"))

if (!file.exists(output_name) | isoverwrite){
  sim <- simdata10(spec,seed=seed_idx[slurm_id])
  max_time <- ceiling(max(sim$data$time))
  if (spec$pl){
    s <- sim$s[sim$s<=max_time]
    k <- length(s)-1
  } else {
    s <- seq(0,max_time,length.out=spec$k+1)
    k <- spec$k
  }
  
  runspec <- list(S=c(spec$bdd_x,spec$bdd_y),maxiter=30,y=sim$data$time,s=s,k=k,
                  delta=sim$delta)

  j <- 1
  runspec <- updateHSGP_v3(runspec,spec,j)
  
  mod_time <- system.time({
    while (!runspec$check[[j]] & j<=runspec$maxiter){
      
      runspec <- runHSGP_v5(runspec,spec,mod,j,seed=seed_idx[slurm_id])
      j <- j+1
      runspec <- updateHSGP_v3(runspec,spec,j)
      
    }
  })
  
  res <- runspec$fit[[j-1]]
  runspec$fit <- NULL
  output <- list(res=res,runspec=runspec,mod_time=mod_time)
  
  saveRDS(output, output_name)
}


