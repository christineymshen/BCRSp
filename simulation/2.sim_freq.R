
specnum <- 1
library(survival)

# folder path to store the simulation results
outfolder <- ""
# folder path which stores "functions.R" and the spec file
infolder <- ""

runnum <- 1
nsim <- 500

source(paste0(infolder,"functions.R"))

if (!file.exists(outfolder)) dir.create(outfolder)
spec <- readRDS(paste0(infolder, "spec_",specnum,".rds"))
n <- spec$n;

set.seed(1)
ncenters <- 9
data_kmeans <- kmeans(spec$d,centers = ncenters)

mod1_niter <- mod2_niter <- numeric(nsim)

for (i in 1:nsim){
  seed <- i
  output_name <- paste0(outfolder,"res",specnum,"_",runnum,"_",seed,".rds")
  
  sim <- simdata10(spec,seed=seed)
  
  if (spec$isModels){
    X <- cbind(spec$X,comorbidity=spec$W)
    p <- spec$p+1
  } else {
    X <- spec$X
    p <- spec$p
  }
  colnames(X)[colnames(X) == "insuranceSelf-Pay"] <- "insuranceSelfPay"
  label <- paste(colnames(X),collapse = "+")
  data1 <- cbind(sim$data,X) %>%
    mutate(event=factor(event))
  
  # maximum number of iteration set to 100 to ensure convergence
  mod1 <- coxph(as.formula(paste0("Surv(time,event) ~",label)), data1, id=c(1:n),
                control=coxph.control(iter.max=100))
  
  data2 <- cbind(data1,group=data_kmeans$cluster) %>% 
    mutate(group=factor(group))
  
  X1 <- model.matrix(~group-1, data=data2)[,-1]
  data2 <- cbind(data2, X1) %>%
    select(-group)
  X <- cbind(X,X1)
  label <- paste(colnames(X),collapse = "+")
  
  # maximum number of iteration set to 100 to ensure convergence
  mod2 <- coxph(as.formula(paste0("Surv(time,event) ~",label)), data2, id=c(1:n),
                control=coxph.control(iter.max=100))
  
  mod1_niter[i] <- mod1$iter
  mod2_niter[i] <- mod2$iter
  
  output <- list(coxph1_m=mod1$coefficients,
                 coxph1_sd=sqrt(diag(mod1$var)),
                 coxph2_m=mod2$coefficients[c(1:p,(p+ncenters):(2*p+ncenters-1))],
                 coxph2_sd=sqrt(diag(mod2$var))[c(1:p,(p+ncenters):(2*p+ncenters-1))])
  
  saveRDS(output, output_name)
  
}

# checking number of iterations used for each run to ensure convergence
range(mod1_niter)
hist(mod1_niter)
which(mod1_niter>=100)
range(mod2_niter)
hist(mod2_niter)
which(mod2_niter>=100)
