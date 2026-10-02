library(tidyverse)
library(rstan)
library(kableExtra)

specnum <- c(2,3,4,9,24,17,18,19)
runnum <- c(1,1,1,1,1,1,1,1)

ns <- length(specnum)
folder <- "application/"
source(paste0(folder,"R/functions.R"))
spec <- readRDS(paste0(folder, "spec/spec",specnum[1],".rds"))

rowspeclabels <- c("Race - Black", "Race - Other", "Sex - Female", "Smoking - Former", 
                   "Smoking - Current", "No partner", "Insurance - Medicaid",
                   "Insurance - Medicare", "Insurance - Other", "Age")
orderedlabels <- c("Age", "Sex - Female", "Race - Black", "Race - Other", "No partner",
                   "Smoking - Former", "Smoking - Current","Insurance - Medicaid",
                   "Insurance - Medicare", "Insurance - Other")
collabels <- c("Low Correlation", "High Correlation", "k=100", "Censoring", "Bumpy",
               "Smooth", "Large tau", "Small tau")

HR1 <- HR2 <- vector("list",length=ns)
for (i in 1:ns){
  fit_i <- readRDS(paste0(folder,"res/spec",specnum[i],"_run",runnum[i],".rds"))$res
  
  p_beta_i <- as.matrix(fit_i,pars="beta")
  HR1[[i]] <- getHR_multi_v3(p_beta_i,rowspeclabels,risktype=1, ordered_labels=orderedlabels)
  HR2[[i]] <- getHR_multi_v3(p_beta_i,rowspeclabels,risktype=2, ordered_labels=orderedlabels)
}

HR1 <- lapply(HR1, function(df) df %>% 
               mutate(HR=paste0(formatC(HR_m,format="f",digits=2), " (",formatC(HR_l,format="f",digits=2),", ", formatC(HR_u,format="f",digits=2), ")")) %>%
               select(-c(label,HR_l,HR_u,HR_m)))

do.call("cbind",HR1) %>%
  `rownames<-`(orderedlabels) %>%
  `colnames<-`(collabels) %>%
  kbl(digits=2, format="latex", row.names = T, align="lcccccccc")

HR2 <- lapply(HR2, function(df) df %>% 
                mutate(HR=paste0(formatC(HR_m,format="f",digits=2), " (",formatC(HR_l,format="f",digits=2),", ", formatC(HR_u,format="f",digits=2), ")")) %>%
                select(-c(label,HR_l,HR_u,HR_m)))

do.call("cbind",HR2) %>%
  `rownames<-`(orderedlabels) %>%
  `colnames<-`(collabels) %>%
  kbl(digits=2, format="latex", row.names = T, align="lcccccccc")
