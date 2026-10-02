# 20260910 In this version, we compare everything to BSp i+s, also include se_diff results

specnum <- c(2,3,4,9,24,17,18,19)
runnum <- c(1,1,1,1,1,1,1,1)

specnum_i <- c(6,7,8,10,25,21,22,23)
runnum_i <- c(1,1,1,1,1,1,1,1)

specnum_np <- c(12,13,14,15,12,12,12,12)
runnum_np <- c(1,1,1,1,1,1,1,1)

folder <- "application/"
library(loo)
library(tidyverse)
library(kableExtra)

ns <- length(specnum)
ncpu <- 6
res_df <- data.frame(num=c(1:2))

orderedlabels <- c("model3", "model2", "model1")

# note that the following code won't always work
# this is based on the fact that I already checked that model 3 always have the lowest WAIC

for (i in 1:ns){
  
  fit <- readRDS(paste0(folder,"res/spec",specnum[i],"_run",runnum[i],".rds"))$res
  fit_i <- readRDS(paste0(folder,"res/spec",specnum_i[i],"_run",runnum_i[i],".rds"))$res
  fit_np <- readRDS(paste0(folder,"res/spec",specnum_np[i],"_run",runnum_np[i],".rds"))$res
  
  loglik <- extract_log_lik(fit,merge_chains = FALSE,parameter_name="loglik")
  r_eff <-   relative_eff(exp(loglik), cores = ncpu)
  loo_is <- loo(loglik, r_eff = r_eff, cores = ncpu)
  
  loglik_i <- extract_log_lik(fit_i,merge_chains = FALSE,parameter_name="loglik")
  r_eff_i <- relative_eff(exp(loglik_i), cores = ncpu)
  loo_i <- loo(loglik_i, r_eff = r_eff_i, cores = ncpu)

  loglik_np <- extract_log_lik(fit_np,merge_chains = FALSE,parameter_name="log_lik")
  r_eff_np <- relative_eff(exp(loglik_np), cores = ncpu)
  loo_np <- loo(loglik_np, r_eff = r_eff_np, cores = ncpu)
  
  waic_res <- loo_compare(waic(loglik_np),waic(loglik_i),waic(loglik))
  
  mod2_res <- paste0(round(waic_res["model2",1]*(-2),1), " (", 
                round(waic_res["model2",2]*2,1), ")")
  mod1_res <- paste0(round(waic_res["model1",1]*(-2),1), " (", 
                     round(waic_res["model1",2]*2,1), ")")
  
  res_df <- cbind(res_df, c(mod2_res,mod1_res))

}

colnames(res_df) <- c("num","low corr","high corr","k=100","censoring",
                      "bumpy","smooth","large tau", "small tau")

# model1: np; model2: i; model3: is
print(res_df %>%
        select(-num) %>%
        knitr::kable() %>%
        kable_styling(full_width = F))


res_df %>%
  select(-num) %>%
  kbl(digits=2, format="latex") 
