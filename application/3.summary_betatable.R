
specnum <- 163

library(tidyverse)
library(kableExtra)

wdir <- ""

folder <- paste0(wdir, "CRS/")
source(paste0(folder,"R/functions.R"))

risklabels <- c("readmission risk", "mortality risk")

# rounding scalar
rs <- 10
# number of intervals
ni <- 4

rowspeclabels <- c("Race - black", "Race - other", "Gender - female", "Smoking - former", "Smoking - current",
                "No partner", "Insurance - comm","Insurance - WCSC", "Insurance - selfpay", "Age",
                "Comorbidity")

  
# read in spec file and results file
spec <- readRDS(paste0(folder, "spec/spec",specnum,".rds"))
fit <- readRDS(paste0(folder,"res/HSGP_spec",specnum,"_run1.rds"))$res

# extract posterior samples
beta <- extract_beta_v2(spec,fit)

## regression coefficients plots

HR <- vector("list", length=spec$m)

for (j in 1:spec$m){
  if (j>1) {
    HR[[j]] <- getHR_multi_v3(beta, rowspeclabels,risktype=j, ordered_labels = HR[[1]]$label)
  } else {
    HR[[j]] <- getHR_multi_v3(beta, rowspeclabels,risktype=j)
  }
}

rownames <- HR[[1]]$label

HR <- lapply(HR, function(df) df %>% 
               mutate(CI=paste0("(",round(HR_l,2),", ", round(HR_u,2), ")")) %>%
               select(-c(label,HR_l,HR_u)))

do.call("cbind",HR) %>%
  `rownames<-`(rownames) %>%
  `colnames<-`(rep(c("Mean", "CI"),2)) %>%
  kbl(digits=2, format="latex", row.names = T, align="lcccc") %>%
  # kbl(digits=2, row.names = T, align="lcccc") %>%
  add_header_above(c(" ", "Readmission Risk"=2, "Mortality Risk"=2))

