slurm_id <- "to be updated"
specnum <- "to be updated"
runnum <- 1

FQrunnum <- 1

isoverwrite <- T
runlabel <- "GP"

folder <- "simulation/"
source(paste0(folder,"R/functions.R"))
outfolder <- paste0(folder,"res/", runlabel)
if (!file.exists(outfolder)) dir.create(outfolder, recursive = T)

spec <- readRDS(paste0(folder, "spec/spec",specnum,".rds"))

stan_file <- "stan/CRS_is_GP3.stan"
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
  
  data <- list(n=spec$n,n_pred=spec$n_pred,p=spec$p,k=k,m=spec$m,
               y=sim$data$time,s=s,delta=sim$delta,X=spec$X,W=spec$W,d=spec$d,
               d_pred=spec$d_pred,a0=spec$a0,b0=spec$b0,a1=spec$a1,b1=spec$b1,
               al=spec$al,bl=spec$bl,tau=spec$tau)
  
  mod_time <- system.time(
    res <- sampling(mod,data=data,chains=4,warmup=spec$nburn,
                    iter=spec$niter,cores=spec$nchain,thin=spec$nthin,
                    seed=seed_idx[slurm_id],include=T,
                    pars=c("beta","beta_w","alpha0","alpha1","l0","l1",
                           "theta_pred","lambda","s2","kappa","theta0","theta1"))
  )
  
  output <- list(res=res,runspec=data,mod_time=mod_time)
  
  saveRDS(output, output_name)
}


