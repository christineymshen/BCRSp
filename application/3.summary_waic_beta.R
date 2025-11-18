
library(loo)
ncpu <- 4
iscv <- T
isbeta <- F
isbeta_plotnp <- F

wdir <- ""

folder <- paste0(wdir, "CRS/")
source(paste0(folder,"R/functions.R"))

outfolder <- ""

# base, low corr, high corr, k=100
specnum_is <- c(163,171,179,187,195)
specnum_i <- specnum_is+4

nrun <- length(specnum_i)
runnums <- rep(1,nrun)

rowspeclabels <- c("Race - black", "Race - other", "Gender - female", "Smoking - former", "Smoking - current",
                   "No partner", "Insurance - comm","Insurance - WCSC", "Insurance - selfpay", "Age",
                   "Comorbidity")

if (is.null(specnum_i)){
  speclabels <- c("with spatial")
} else {
  speclabels <- c("spatial i+s", "spatial i")
}

if (isbeta_plotnp) speclabels <- c(speclabels, "no spatial")

for (i in 1:nrun){
  
  spec_label <- paste0(specnum_is[i],"_",specnum_i[i])
  
  spec_i <- readRDS(paste0(folder, "spec/spec",specnum_i[i],".rds"))
  spec_is <- readRDS(paste0(folder, "spec/spec",specnum_is[i],".rds"))
  
  # read in spec file and results file
  spec_is <- readRDS(paste0(folder, "spec/spec",specnum_is[i],".rds"))
  fit_is <- readRDS(paste0(folder,"res/HSGP_spec",specnum_is[i],"_run",runnums[i],".rds"))$res
  
  spec_i <- readRDS(paste0(folder, "spec/spec",specnum_i[i],".rds"))
  fit_i <- readRDS(paste0(folder,"res/HSGP_spec",specnum_i[i],"_run",runnums[i],".rds"))$res

  if (file.exists(paste0(folder,"res/CRS_np_spec", specnum_i[i], "_run1.rds"))){
    specnum_np <- specnum_i[i]
  } else {
    specnum_np <- specnum_is[i]
  }
  
  fit_np <- readRDS(paste0(folder,"res/CRS_np_spec",specnum_np,"_run",runnums[i],".rds"))$res

  outcome <- spec_i$outcome
  
  ### CV results
  if (iscv){
    
    if (!file.exists(paste0(folder,"res/r5/waic_spec",spec_label,".rds"))){
      loglik_np <- extract_log_lik(fit_np,merge_chains = FALSE,parameter_name="log_lik")
      loglik_is <- extract_log_lik(fit_is,merge_chains = FALSE,parameter_name="loglik")
      r_eff_np <- relative_eff(exp(loglik_np), cores = ncpu)
      r_eff_is <- relative_eff(exp(loglik_is), cores = ncpu)
      loo_np <- loo(loglik_np, r_eff = r_eff_np, cores = ncpu)
      loo_is <- loo(loglik_is, r_eff = r_eff_is, cores = ncpu)
      
      if (!is.null(specnum_i)) {
        loglik_i <- extract_log_lik(fit_i,merge_chains = FALSE,parameter_name="loglik")
        r_eff_i <- relative_eff(exp(loglik_i), cores = ncpu)
        loo_i <- loo(loglik_i, r_eff = r_eff_i, cores = ncpu)
      } 
      
      
      if (is.null(specnum_i)){
        waic_res <- loo_compare(waic(loglik_np),waic(loglik_is))
        ploocv_res <- loo_compare(loo_np,loo_is)
        outputlabel <- specnum_is[i]
      } else {
        waic_res <- loo_compare(waic(loglik_np),waic(loglik_i),waic(loglik_is))
        ploocv_res <- loo_compare(loo_np,loo_i,loo_is)
        outputlabel <- paste0(specnum_is[i],"_",specnum_i[i])
        
      }
      
      saveRDS(waic_res, paste0(folder,"res/r5/waic_spec",outputlabel,".rds"))
      saveRDS(ploocv_res, paste0(folder,"res/r5/ploocv_spec",outputlabel,".rds"))

    }
    
    res <- data.frame(num=c(1:3))
    
    waic_res <- readRDS(paste0(folder,"res/r5/waic_spec",spec_label,".rds"))
    ploocv_res <- readRDS(paste0(folder,"res/r5/ploocv_spec",spec_label,".rds"))
    
    res <- cbind(res, paste0(rownames(waic_res),", ", round(waic_res[,7],2)))
    res <- cbind(res, paste0(rownames(waic_res),", ", round(waic_res[,5],2)))
    res <- cbind(res, paste0(rownames(ploocv_res),", ", round(ploocv_res[,3],2)))
    
    colnames(res) <- c("num","waic","p_waic","pareto-smoothed loocv")
    
    # model1: np; model2: i; model3: is
    print(res %>%
      select(-num) %>%
      knitr::kable() %>%
      kable_styling(full_width = F))
  }
  
  if (isbeta){
    ### beta figures
    # extract posterior samples
    
    beta_np <- as.matrix(fit_np,pars="beta")
    beta_is <- extract_beta_v2(spec_is,fit_is)
    if (!is.null(specnum_i)) beta_i <- extract_beta_v2(spec_i,fit_i)
    
    ## regression coefficients plots
    png(file = paste0(outfolder,"fig/spec", specnum_is[i], ifelse(is.null(specnum_i), "", paste0("_",specnum_i[i])), "_beta.png"), width=800, height=600)
    par(mfrow=c(1,2),mar=c(3,1,2,1), mgp=c(1.75, 0.75, 0), oma=c(1,9,1,1))
    for (j in 1:m){
      HR_is <- getHR_multi_v3(beta_is, rowspeclabels,risktype=j)
      values <- HR_is[,3]
      dfs <- list(HR_is)
      
      if (isbeta_plotnp){
        HR_np <- getHR_multi_v3(beta_np, rowspeclabels,risktype=j)
        values <- c(values, HR_np[,3])
        dfs <- append(dfs,list(HR_np))
      }
      
      if (!is.null(specnum_i)) {
        HR_i <- getHR_multi_v3(beta_i, rowspeclabels,risktype=j)
        values <- c(values, HR_i[,3])
        dfs <- append(dfs,list(HR_i))
      }
      
      xmax <- ceiling(max(values)*rs)/rs
      
      if (j==1){
        plotHRs_v3(dfs, risklabels[j], cex=1.2, xlim=c(0,xmax), 
                   legend_labels = speclabels, legend_position = "topleft")
        label_record <- HR_is$label
      } else {
        plotHRs_v3(dfs, risklabels[j], label_record, yaxt=F, cex=1.2, xlim=c(0,xmax))
      }
      
    }
    dev.off() 
    
    
  }
}

