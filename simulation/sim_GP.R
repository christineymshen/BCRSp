# slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# specnum <- as.integer(Sys.getenv("specnum"))

slurm_id <- 1
specnum <- 1

outfolder <- ""
infolder <- ""

library(rstan)

runnum <- 1 # in case multiple runs were needed for the same spec
nsim <- 2 # number of simulations per slurm_id
isoverwrite <- T # whether overwrite existing saved run results
source(paste0(infolder,"functions.R"))

if (!file.exists(outfolder)) dir.create(outfolder)

spec <- readRDS(paste0(infolder, "spec_",specnum,".rds"))
n <- spec$n
start_seed <- (slurm_id-1)*nsim
mod <- stan_model(paste0(infolder, "CRS_is_GP2.stan"))

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
    
    data <- list(n=n,n_pred=spec$n_pred,p=spec$p,k=k,m=spec$m,
                 y=sim$data$time,s=s,delta=sim$delta,X=spec$X,W=spec$W,d=spec$d,
                 d_pred=spec$d_pred,kappa0=spec$kappa0,a=spec$a,b=spec$b)
    
    mod_time <- system.time(
      res <- sampling(mod,data=data,chains=spec$nchain,warmup=spec$nburn,
                      iter=spec$niter,cores=spec$nchain,thin=spec$nthin,
                      seed=seed,refresh=0,include=T,
                      pars=c("beta","beta_w","alpha0","alpha1","l0","l1",
                             "theta_pred","lambda","s2","kappa1"))
    )
    
    output <- list(res=res,runspec=data,mod_time=mod_time)
    
    saveRDS(output, output_name)
  }
}
