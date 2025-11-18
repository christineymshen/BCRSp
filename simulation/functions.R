library(tidyverse)
library(matrixStats)
library(doMC)
library(viridis)
library(ggpubr)

# To compute Matern covariance matrix
matern_cov <- function(d,alpha,nu,l){
  if (d==0){
    alpha^2
  } else {
    x <- sqrt(2*nu)*d/l
    alpha^2*2^(1-nu)/gamma(nu)*x^nu*besselK(x,nu)
  }
}
vec_matern_cov <- Vectorize(matern_cov)

# This function turns vectors computed from the dist function back to a matrix
distvec_to_mat <- function(d_vec,n){
  
  d_mat <- matrix(0,n,n)
  
  d_mat[lower.tri(d_mat)] <- d_vec
  d_mat <- d_mat + t(d_mat)
  
  diag(d_mat) <- 1
  return(d_mat)
}

# Simulate competing risks datasets
simdata10 <- function(spec,seed=1){
  
  # pc: proportion of censored data points, default at 0.2
  # pl: piecewise linear, TURE or FALSE
  
  n <- spec$n; m <- spec$m; nc <- floor(spec$pc*n) # number censored
  pl <- spec$pl
  
  if (spec$isModeli & spec$isModels){
    theta0_true <- cbind(spec$f_true[1,1,],spec$f_true[2,1,])[1:n,]
    theta1_true <- cbind(spec$f_true[1,2,],spec$f_true[2,2,])[1:n,]
    Xbeta <- spec$X %*% spec$beta + theta0_true + (theta1_true + outer(rep(1,n), spec$beta_w)) * spec$W
  } else if (spec$isModeli & !spec$isModels){
    theta0_true <- cbind(f_true[1,1,],f_true[2,1,])[1:n,]
    Xbeta <- spec$X %*% spec$beta + theta0_true
  } else if (!spec$isModeli & spec$isModels) {
    theta1_true <- cbind(f_true[1,1,],f_true[2,1,])[1:n,]
    Xbeta <- spec$X %*% spec$beta + (theta1_true + outer(rep(1,n), spec$beta_w)) * spec$W
  }
  # deliberately put "shift" together with Xbeta, for numerical stability
  Xbeta <- Xbeta + outer(rep(1,n),spec$shift) 
  
  set.seed(seed)
  Us <- matrix(runif(n*m),n,m)
  U_Xbeta <- -log(Us)/exp(Xbeta)
  
  lambda_f <- partial(spec$lambda_f_sim,spec=spec)
  vec_lambda_f <- Vectorize(lambda_f)
  if (pl){
    s <- spec$s
    k <- length(s)-1
    lambda_s <- matrix(0,k+1,m)
    # effectively, the rate at s[2] is the rate applied to the first interval
    for (j in 1:m){
      lambda_s[1:k,j] <- vec_lambda_f(tail(s,k),j)
      lambda_s[k+1,j] <- lambda_s[k,j]
    }
    H_s <- matrix(0,k+1,m)
    H_s[2:(k+1),] <- lambda_s[1:k,] * (tail(s,k)-head(s,k))
    H_s <- apply(H_s,2,cumsum)
    
    idx <- vapply(c(1:m),function(j) findInterval(U_Xbeta[,j],H_s[,j]), integer(n))
    T_latent <- vapply(c(1:m),function(j) (U_Xbeta[,j]-H_s[idx[,j],j])/lambda_s[idx[,j],j]
                       + s[idx[,j]], numeric(n))
  } else {
    invert_cdf <- partial(spec$invert_cdf,spec=spec)
    vec_invert_cdf <- Vectorize(invert_cdf)
    T_latent <- matrix(0,n,m)
    for (j in 1:m){
      T_latent[,j] <- vec_invert_cdf(U_Xbeta[,j],j)
    }
  }
  
  # get censored data
  idx_c <- sample(1:n, nc)
  
  ctime_upper <- rowMins(T_latent[idx_c,])
  ctime <- runif(nc,0,ctime_upper)
  
  time <- rowMins(T_latent)
  time[idx_c] <- ctime
  event <- apply(T_latent,1,which.min)
  event[idx_c] <- 0
  data <- data.frame(cbind(time,event))
  
  delta <- matrix(0,n,m)
  idx <- cbind(c(1:n), event)[event!=0,]
  delta[idx] <- 1
  
  res <- list(data=data,delta=delta)
  if (pl) res$s <- s
  
  return(res)
}

### HSGP helper functions

# these are functions to be used to tune hyperparameters for HSGP under matern 32 kernel
m_QE <- function(c,l,S) ceiling(3.42 * c / (l/S))
l_QE <- function(c,m,S) round(S * 3.42 * c / m, 3) 

c_vs_l_QE <- function(l,S){
  c =  4.5*l/S
  if(c < 1.2)
    c = 1.2
  c
}

diagnostic <- function(l,l_hat) l_hat + 0.01 >= l


updateHSGP_v2 <- function(rs,ss,i){
  S <- rs$S
  d <- 2
  
  if (rs$isModeli & rs$isModels & !rs$isJointwithin) {
    nGP <- 2
    l_labels <- c("l0","l1")
  } else {
    nGP <- 1
    if (rs$isJointwithin) {
      l_labels <- "l"
    } else if (rs$isModeli) {
      l_labels <- "l0"
    } else if (rs$isModels) {
      l_labels <- "l1"
    }
  }
  
  m <- sum(rs$isModelrt)
  
  # we assume each GP has only one lengthscale parameter, for all the dimensions
  # but HSGP setup would require one lengthscale for each dimension during configuration
  # because the length of each dimension could be different
  
  ci <- mi <- lij <- array(0,dim=c(nGP,d,m))
  li <- li_hat <- matrix(0,nGP,m)
  diagi <- check <- matrix(F,nGP,m)
  
  if (i==1){
    li <- matrix(0.5,nGP,m)
    for (k in 1:nGP){
      for (j in 1:d){
        for (l in 1:m){
          ci[k,j,l] <- c_vs_l_QE(l=li[k,l], S=S[j])
          mi[k,j,l] <- m_QE(c=ci[k,j,l], l=li[k,l], S=S[j])          
        }
      }
    }   
  } else {
    for (k in 1:nGP){
      li_hat[k,] <- round(summary(rs$fit[[i-1]], 
                                  pars = l_labels[k], probs = c(0.05, 0.95))$summary[,1], 2)
      diagi[k,] <- diagnostic(rs$l[[i-1]][k,], li_hat[k,])
    }
    
    if (i==2){
      
      for (k in 1:nGP){
        for (l in 1:m){
          if (diagi[k,l]){
            mi[k,,l] <- rs$m[[i-1]][k,,l] + 2
            for (j in 1:d){
              ci[k,j,l] <- c_vs_l_QE(l=li_hat[k,l], S=S[j])
              lij[k,j,l] <- l_QE(ci[k,j,l], m=mi[k,j,l], S=S[j])
            }
            li[k,l] <- max(lij[k,,l])
          } else {
            li[k,l] <- li_hat[k,l]
            for (j in 1:d){
              ci[k,j,l] <- c_vs_l_QE(l=li[k,l], S=S[j])
              mi[k,j,l] <- m_QE(c=ci[k,j,l], l=li[k,l], S=S[j])
            }
          }
        }
      }
      
    } else {
      
      for (k in 1:nGP){
        for (l in 1:m){
          if (diagi[k,l] & !rs$diag[[i-2]][k,l]){
            mi[k,,l] <- rs$m[[i-1]][k,,l] + 2
            for (j in 1:d){
              ci[k,j,l] <- c_vs_l_QE(l=li_hat[k,l], S=S[j])
              lij[k,j,l] <- l_QE(ci[k,j,l], m=mi[k,j,l], S=S[j])
            }
            li[k,l] <- max(lij[k,,l])
          } else if (diagi[k,l] & rs$diag[[i-2]][k,l]){
            check[k,l] <- T
            li[k,l] <- rs$l[[i-1]][k,l]
            ci[k,,l] <- rs$c[[i-1]][k,,l]
            mi[k,,l] <- rs$m[[i-1]][k,,l]
          } else if(!diagi[k,l]){
            li[k,l] <- li_hat[k,l]
            for (j in 1:d){
              ci[k,j,l] <- c_vs_l_QE(l=li[k,l], S=S[j])
              # for the first 5 runs before things stabilize, we don't apply the flooring
              if (i>5){
                mi[k,j,l] <- max(m_QE(c=ci[k,j,l], l=li[k,l], S=S[j]), rs$m[[i-1]][k,j,l])
              } else {
                mi[k,j,l] <- m_QE(c=ci[k,j,l], l=li[k,l], S=S[j])
              }
            }       
          }
        }
      }
    }
    
    rs$l_hat[[i-1]] <- li_hat
    rs$diag[[i-1]] <- diagi
  }
  
  if (i==1){
    rs$fit <- list()
    rs$l <- list() # lengthscale
    rs$c <- list()  # boundary factor
    rs$m <- list()  # number of basis functions
    rs$l_hat <- list()  # HSGP lengthscale estimate
    rs$diag <- list()
    rs$check <- list()
  }
  rs$check[[i]] <- all(check)
  rs$l[[i]] <- li
  rs$c[[i]] <- ci
  rs$m[[i]] <- mi
  
  return(rs)
}

runHSGP_v3 <- function(rs,ss,mod,i,seed=1,testshort=F,testfull=F){
  
  if (testshort){
    t_iter <- 200; t_burn <- 100
    refresh <- max(10, t_iter/10)
  } else {
    t_iter <- ss$niter; t_burn <- ss$nburn
    
    if (testfull){
      refresh <- max(10, t_iter/10)
    } else {
      refresh <- 0
    }
  }
  
  if (sum(rs$isModelrt)==2 & rs$isJointwithin){
    data <- list(n=ss$n,n_pred=ss$n_pred,p=ss$p,k=rs$k,y=rs$y,s=rs$s,
                 delta=rs$delta,X=ss$X,W=ss$W,d=ss$d,d_pred=ss$d_pred,a=ss$a,
                 b=ss$b,kappa0=ss$kappa0,L=rs$c[[i]][1,,]*rs$S,M=rs$m[[i]][1,,])
    
    rs$fit[[i]] <- sampling(mod, data = data, init=0.5, chains = ss$nchain, 
                            warmup = t_burn,iter = t_iter, core=ss$nchain, 
                            thin = ss$nthin, seed = seed, include=T, refresh=refresh,
                            pars = c("alpha0","l","alpha1","rho","beta","beta_w",
                                     "s2","kappa1","lambda","theta_pred","alpha0_adj",
                                     "alpha1_adj")) 
    
  } else if (rs$isModels & !rs$isModeli){
    data <- list(n=ss$n,n_pred=ss$n_pred,p=ss$p,k=rs$k,m=ss$m,y=rs$y,s=rs$s,
                 delta=rs$delta,X=ss$X,W=ss$W,d=ss$d,d_pred=ss$d_pred,a=ss$a,
                 b=ss$b,kappa0=ss$kappa0,L=rs$c[[i]][1,,]*rs$S,M=rs$m[[i]][1,,])
    
    rs$fit[[i]] <- sampling(mod, data = data, init=0.5, chains = ss$nchain, 
                            warmup = t_burn,iter = t_iter, core=ss$nchain, 
                            thin = ss$nthin, seed = seed, include=T, refresh=refresh,
                            pars = c("l1","alpha1","beta","beta_w","s2","kappa1",
                                     "lambda","theta_pred","alpha1_adj"))  
    
  } else if (rs$isModeli & rs$isModels){
    data <- list(n=ss$n,n_pred=ss$n_pred,p=ss$p,k=rs$k,y=rs$y,s=rs$s,
                 delta=rs$delta,X=ss$X,W=ss$W,d=ss$d,d_pred=ss$d_pred,a=ss$a,
                 b=ss$b,kappa0=ss$kappa0,L0=rs$c[[i]][1,,]*rs$S,M0=rs$m[[i]][1,,],
                 L1=rs$c[[i]][2,,]*rs$S,M1=rs$m[[i]][2,,])
    
    rs$fit[[i]] <- sampling(mod, data = data, init=0.5, chains = ss$nchain, 
                            warmup = t_burn,iter = t_iter, core=ss$nchain, 
                            thin = ss$nthin, seed = seed, include=T, refresh=refresh,
                            pars = c("l0","alpha0","l1","alpha1","beta","beta_w",
                                     "s2","kappa1","lambda","theta_pred",
                                     "alpha0_adj","alpha1_adj"))      
    
  }
  
  return(rs)
}

# Get summary statistics across results on all simulate datasets for GP / HSGP
get_stats_v5 <- function(infolder,resfolder,outfolder,runtype,specnum,runnum,nsim,n_cpus=6){
  
  registerDoMC(cores = n_cpus)
  isHSGP <- runtype=="HSGP"
  
  spec <- readRDS(paste0(infolder, "spec_",specnum,".rds"))

  isModeli <- spec$isModeli
  isModels <- spec$isModels
  isJointwithin <- spec$isJointwithin
  isModelrt <- spec$isModelrt
  run_s_idx <- spec$s_idx
  
  path <- paste0(resfolder, "res",specnum,"_",runnum,"_")
  
  n <- spec$n; m <- spec$m; n_pred <- spec$n_pred
  f_true_check <- spec$f_true_check # m x nGP x n_pred
  
  nGP <- isModeli+isModels
  
  if (isModeli) {
    f0_true <- t(f_true_check[,1,])
  }
  if (isModels){
    W <- spec$W
    f1_true <- t(f_true_check[,run_s_idx,])
  } 
  
  para_res <- get_labels_v2(spec)
  t_para <- para_res$para_val
  pars <- para_res$para_name
  
  n_para <- length(t_para)
  n_psim <- (spec$niter-spec$nburn)/spec$nthin*spec$nchain
  n_psim1 <- rep(1,n_psim)
  
  lambda_f <- partial(spec$lambda_f_check,spec=spec)
  vec_lambda_f <- Vectorize(lambda_f)
  
  if (spec$pl){
    s_all <- spec$s
    k_all <- length(s_all)-1
    t_lambda_all <- matrix(0,k_all,m)
    for (j in 1:m){
      t_lambda_all[,j] <- vec_lambda_f(tail(s_all,k_all),j)
    }
  } else {
    k <- spec$k
    # for storage
    t_lambda <- matrix(0,k,m)
  }
  
  res <- foreach(i=1:nsim) %dopar% {
    
    res_i <- readRDS(paste0(path,i,".rds"))
    fit_i <- res_i$res
    s <- res_i$runspec$s
    
    # obtain true lambda's, need to be put inside the for loop because of s
    if (spec$pl){
      k <- length(s)-1
      t_lambda <- t_lambda_all[1:k,]
    } else {
      for (j in 1:m){
        for (t in 1:k){
          t_lambda[t,j] <- integrate(vec_lambda_f,s[t],s[t+1],j=j)$value/(s[t+1]-s[t])
        }
      }
    }
    
    # extract posterior samples for parameters other than theta and lambda
    para <- as.matrix(fit_i,pars=pars)
    
    # n_psim by m by n_pred by nGP
    theta_pred <- array(rstan::extract(fit_i,"theta_pred")$theta_pred,
                        dim=c(n_psim,m,n_pred,nGP))
    
    if (isModeli){
      theta0_pred <- theta_pred[,,,1]
      
      theta0_pred_level <- xbeta_level <- matrix(0,n_psim,m)
      theta0_sd <- theta0_std_diff <- theta0_m <- theta0_coverage <- theta0_rmse <- matrix(0,n_pred,m)
      
      for (j in 1:m){
        theta0_pred_level[,j] <- rowMeans(theta0_pred[,j,])
        xbeta_level[,j] <- colMeans(tcrossprod(spec$X,para[,((j-1)*spec$p+1):(j*spec$p)]))
        theta0_pred[,j,] <- theta0_pred[,j,] - theta0_pred_level[,j]
        theta0_sd[,j] <- apply(theta0_pred[,j,],2,sd)
        theta0_m[,j] <- colMeans(theta0_pred[,j,])
        theta0_m_diff_j <- theta0_m[,j]-f0_true[,j]
        theta0_std_diff[,j] <- theta0_m_diff_j / theta0_sd[,j]
        theta0_CI_j <- t(apply(theta0_pred[,j,],2,quantile,c(0.025,0.975)))
        theta0_coverage[,j] <- f0_true[,j] > theta0_CI_j[,1] & f0_true[,j] < theta0_CI_j[,2]
        theta0_rmse[,j] <- sqrt(colMeans(theta0_pred[,j,] - outer(n_psim1,f0_true[,j]))^2)
        
      }
    }
    
    if (isModels){
      theta1_pred <- theta_pred[,,,run_s_idx]
      
      theta1_pred_level <- matrix(0,n_psim,m)
      p_beta_w <- as.matrix(fit_i,pars="beta_w")
      theta1_pred <- vapply(c(1:n_pred), function(i) theta1_pred[,,i]+p_beta_w, matrix(0,n_psim,m)) #npsim x m x n_pred
      
      theta1_sd <- theta1_m <- theta1_std_diff <- theta1_coverage <- theta1_rmse <- matrix(0,n_pred,m)
      
      for (j in 1:m){
        theta1_pred_level[,j] <- rowMeans(theta1_pred[,j,])
        theta1_sd[,j] <- apply(theta1_pred[,j,],2,sd)
        theta1_m[,j] <- colMeans(theta1_pred[,j,])
        theta1_m_diff_j <- theta1_m[,j]-f1_true[,j]
        theta1_std_diff[,j] <- theta1_m_diff_j / theta1_sd[,j]
        theta1_CI_j <- t(apply(theta1_pred[,j,],2,quantile,c(0.025,0.975)))
        theta1_coverage[,j] <- f1_true[,j] > theta1_CI_j[,1] & f1_true[,j] < theta1_CI_j[,2]
        theta1_rmse[,j] <- sqrt(colMeans(theta1_pred[,j,] - outer(n_psim1,f1_true[,j]))^2)
      }
      
      idx <- which(str_detect(colnames(para),"beta_w"))
      para[,idx] <- theta1_pred_level
    }
    
    lambda <- rstan::extract(fit_i, pars="lambda")$lambda
    lambda_m <- lambda_std_diff <- lambda_rmse <- lambda_coverage <- matrix(0,k,m)
    
    for (j in 1:m){
      if (isModeli)
        lambda[,j,] <- lambda[,j,] * exp(theta0_pred_level[,j]+xbeta_level[,j])
      
      lambda_m[,j] <- colMeans(lambda[,j,])
      lambda_sd_j <- apply(lambda[,j,],2,sd)
      lambda_std_diff[,j] <- (lambda_m[,j] - t_lambda[,j]) / lambda_sd_j
      lambda_rmse[,j] <- sqrt(colMeans( (lambda[,j,] - outer(n_psim1,t_lambda[,j]))^2 ))
      lambda_CI_j <- t(apply(lambda[,j,],2,quantile,c(0.025,0.975)))
      lambda_coverage[,j] <- t_lambda[,j] > lambda_CI_j[,1] & t_lambda[,j] < lambda_CI_j[,2]
    }
    
    para_m <- colMeans(para)
    para_sd <- apply(para,2,sd)
    para_std_diff <- (para_m-t_para)/para_sd
    para_m_diff <- para_m-t_para
    para_rmse <- sqrt(colMeans( (para - outer(n_psim1,t_para))^2 ))
    para_CI <- t(apply(para,2,quantile,c(0.025,0.975)))
    para_coverage <- t_para>para_CI[,1] & t_para <para_CI[,2]
    
    # store ESS
    ESS <- summary(fit_i)$summary[,"n_eff"][colnames(para)] # in this way, the order of the parameters are aligned
    
    # store runtime
    RT <- res_i$mod_time[3]
    
    if (isHSGP) {
      
      alpha_adj <- NULL
      idx <- numeric(0)
      
      if (isModeli){
        alpha_adj <- cbind(alpha_adj, as.matrix(fit_i,pars="alpha0_adj"))
        idx <- c(idx,which(str_detect(colnames(para),"alpha0")))
      }
      
      if (isModels){
        alpha_adj <- cbind(alpha_adj,as.matrix(fit_i,pars="alpha1_adj"))
        idx <- c(idx,which(str_detect(colnames(para),"alpha1")))
      }
      
      alpha_adj_m <- colMeans(alpha_adj)
      alpha_adj_sd <- apply(alpha_adj,2,sd)
      alpha_adj_std_diff <- (alpha_adj_m-t_para[idx])/alpha_adj_sd
      alpha_adj_rmse <- sqrt(colMeans( (alpha_adj - outer(n_psim1,t_para[idx]))^2 ))
      alpha_adj_CI <- t(apply(alpha_adj,2,quantile,c(0.025,0.975)))
      alpha_adj_coverage <- t_para[idx]>alpha_adj_CI[,1] & t_para[idx]<alpha_adj_CI[,2]
      
    }
    
    res_i_tmp <- list(ESS=ESS,RT=RT,para_std_diff=para_std_diff,para_rmse=para_rmse,
                      para_coverage=para_coverage,para_m_diff=para_m_diff,
                      lambda_std_diff=lambda_std_diff,lambda_m=lambda_m,
                      lambda_rmse=lambda_rmse,lambda_coverage=lambda_coverage)
    
    if (isModeli)
      res_i_tmp <- append(res_i_tmp, list(theta0_m=theta0_m,theta0_rmse=theta0_rmse,
                                          theta0_std_diff=theta0_std_diff,theta0_sd=theta0_sd,
                                          theta0_coverage=theta0_coverage))
    if (isModels)
      res_i_tmp <- append(res_i_tmp, list(theta1_m=theta1_m,theta1_rmse=theta1_rmse,
                                          theta1_std_diff=theta1_std_diff,theta1_sd=theta1_sd,
                                          theta1_coverage=theta1_coverage))
    
    if (isHSGP)
      res_i_tmp <- append(res_i_tmp, list(alpha_adj_rmse=alpha_adj_rmse,
                                          alpha_adj_std_diff=alpha_adj_std_diff,
                                          alpha_adj_coverage=alpha_adj_coverage))
    
    if (!spec$pl)
      res_i_tmp$s=s
    
    return(res_i_tmp)
  }
  
  return(res)
}

# Get summary statistics across results on all simulate datasets for frequentist runs
get_stats_freq_v2 <- function(infolder,resfolder,outfolder,specnum,runnum,nsim,n_cpus=6){
  
  registerDoMC(cores = n_cpus)
  spec <- readRDS(paste0(infolder, "spec_",specnum,".rds"))
  path <- path <- paste0(resfolder, "res",specnum,"_",runnum,"_")
  
  n <- spec$n; m <- spec$m; n_pred <- spec$n_pred
  
  p <- spec$p + 1
  para_res <- get_labels_v2(spec)
  t_para <- c(para_res$para_val[1:(p-1)], 
              para_res$para_val[2*(p-1)+1],
              para_res$para_val[p:(2*(p-1))],
              para_res$para_val[2*p])
  
  res <- foreach(i=1:nsim) %dopar% {
    res_i <- readRDS(paste0(path,i,".rds"))
    
    coxph1_std_diff <- (res_i$coxph1_m-t_para)/res_i$coxph1_sd
    coxph2_std_diff <- (res_i$coxph2_m-t_para)/res_i$coxph2_sd
    coxph1_m_diff <- res_i$coxph1_m-t_para
    coxph2_m_diff <- res_i$coxph2_m-t_para
    
    summary_i <- list(coxph1_std_diff=coxph1_std_diff,
                      coxph2_std_diff=coxph2_std_diff,
                      coxph1_m_diff=coxph1_m_diff,
                      coxph2_m_diff=coxph2_m_diff)
    return(summary_i)
  }
  
  return(res)
}

# helper function to get parameter labels
get_labels_v2 <- function(spec){
  
  m <- spec$m
  
  isModeli <- spec$isModeli
  isModels <- spec$isModels
  isJointwithin <- spec$isJointwithin
  
  ### labels and true values for parameters other than theta and lambda
  para_val <- c(spec$beta)
  para_name <- "beta"
  para_latex <- sprintf("b[%d%d]",rep(1:spec$p,m),rep(1:m,each=spec$p))
  
  # these values are always simulated separately as beta_w
  # but they'll be modeled under "beta" if we don't model spatial slopes
  para_val <- c(para_val,spec$beta_w)
  if (isModels){
    para_name <- c(para_name,"beta_w")
    para_latex <- c(para_latex, sprintf("b[w%d]",c(1:m)))
  }
  
  if (isModeli){
    para_val <- c(para_val,spec$alpha[,1])
    para_name <- c(para_name,"alpha0")
    para_latex <- c(para_latex,sprintf("alpha[0*%d]",c(1:m)))
  }
  if (isModels){
    para_val <- c(para_val,spec$alpha[,spec$s_idx])
    para_name <- c(para_name,"alpha1")
    para_latex <- c(para_latex,sprintf("alpha[1%d]",c(1:m)))
  }  
  
  if (!spec$isLMC){
    para_val <- c(para_val,spec$l)
    if (isJointwithin){
      para_name <- c(para_name,"l")
      para_latex <- c(para_latex,sprintf("l[%d]",c(1:m)))
    } else {
      if (isModeli){
        para_name <- c(para_name,"l0")
        para_latex <- c(para_latex,sprintf("l[0*%d]",c(1:m)))
      }
      if (isModels){
        para_name <- c(para_name,"l1")
        para_latex <- c(para_latex,sprintf("l[1%d]",c(1:m)))
      }
    }
    
    if (isJointwithin){
      para_val <- c(para_val,spec$rho)
      para_name <- c(para_name,"rho")
      para_latex <- c(para_latex,sprintf("rho[%d]",c(1:m)))
    }
  }
  
  list(para_val=para_val,para_name=para_name,para_latex=para_latex)
}

# helper functions for spatial surface plots
to_get_plot_v5 <- function(plotf,title,subtitles,rt,dfs,legend=NULL,
                           myscale=NULL,inserts=NULL){
  
  ndf <- length(dfs)
  
  legend <- ifelse(is.null(legend),"",legend)
  
  j <- rt
  
  if (is.null(myscale)){
    m_scale <- get_plot_scale_v3(rs,ni,j,dfs)
  } else {
    m_scale <- myscale
  }
  
  ps <- vector("list", ndf)
  
  for (i in 1:ndf){
    ps[[i]] <- plotf(dfs[[i]][,j], subtitles[i],legend,m_scale)
  }
  
  plot <- ggarrange(plotlist=ps,ncol=ndf,common.legend=T,legend="bottom")
  
  if (!is.null(inserts)){
    plot <- ggarrange(inserts[[j]], plot, ncol=2, widths=c(1,2))
  }
  if (title!="")
    plot <- annotate_figure(plot, fig.lab = title, fig.lab.size = 16,
                            fig.lab.pos="top.left", fig.lab.face="bold")   
  
  print(plot)
  
}

# plot spatial surfaces
get_plot1_v6 <- function(df,title,fill,m_scale,geo_obs=NULL,center=0){
  
  if (!is.null(geo_obs))
    n <- dim(geo_obs)[1]
  
  ncolor <- 10
  if (m_scale$ll<center & m_scale$uu>center){
    if (center-m_scale$ll < m_scale$uu-center) {
      ncolor1 <- max(round(((center-m_scale$ll)/(m_scale$uu-center))*ncolor,0),1)
      ncolor2 <- ncolor
    } else {
      ncolor1 <- ncolor
      ncolor2 <- max(round(((m_scale$uu-center)/(center-m_scale$ll))*ncolor,0),1)
    }
    
    color_scale <- c(tail(viridis(15, option = "G"),ncolor1), "white", 
                     rev(tail(viridis(15, option = "A"),ncolor2)))
    color_values <- scales::rescale(seq(m_scale$ll,m_scale$uu, length.out=length(color_scale)))
    if (!is.null(geo_obs))
      grey_idx <- rep(T,n)
    
  } else {
    color_values <- c(0,1)
    if (m_scale$ll>=center){
      color_scale <- rev(tail(viridis(15, option = "A"),ncolor))
    } else if (m_scale$uu<=center){
      color_scale <- tail(viridis(15, option = "G"),ncolor)
    }
    if (!is.null(geo_obs))
      grey_idx <- rep(F,n)
  }
  
  p <- ggplot(get_df_v2(df)) +
    labs(title=title,x="",y="",fill=fill)+
    geom_tile(aes(x=x,y=y,fill=value)) +
    scale_fill_gradientn(colors=color_scale, values = color_values,
                         breaks=round(seq(m_scale$bl,m_scale$bu,by=m_scale$bby),2),
                         limits=c(m_scale$ll,m_scale$uu)) + 
    theme_bw() + 
    theme(panel.border = element_blank(),
          legend.position="bottom",
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(margin = margin(b = -5)),
          plot.margin = margin(t=15),   # overall plot margin
          legend.margin = margin(t = -15, unit = "pt"))
  
  if (is.null(geo_obs)){
    p
  } else {
    p + geom_point(data=geo_obs[grey_idx,],aes(x=x,y=y),col="grey20",size=0.4,alpha=0.5) + 
      geom_point(data=geo_obs[!grey_idx,],aes(x=x,y=y),col="white",size=0.4)
  }
  
}

# helper function to get boxplots for baseline hazard rates
to_get_plot3_v2 <- function(dfs,axis_label,legend_label,title,
                            ylim=NULL,xbreaks=NULL, ylab=NULL){
  
  if (is.null(ylim)){
    m_scale <- get_plot_scale2(dfs,rs)
    if (m_scale$ll>=0){
      ylim <- c(0,m_scale$ul)
    } else {
      ylim <- c(-m_scale$abs_ul,m_scale$abs_ul)
    }
  }
  
  get_plot3_v2(dfs,axis_label,legend_label,title,ylim,xbreaks,ylab)
}

# boxplot for baseline hazard rate results
get_plot3_v2 <- function(dfs_raw,paranames,labels,title,
                         ylim=NULL,xbreaks=NULL,ylab=NULL){
  
  npara <- length(paranames)
  ndf <- length(dfs_raw)
  dfs <- vector("list",length=ndf)
  if (is.null(ylab)) ylab <- ""
  
  if (is.null(xbreaks)) xbreaks <- c(1:length(paranames))
  
  for (i in 1:ndf){
    dfs[[i]] <- data.frame(dfs_raw[[i]]) %>%
      `colnames<-`(paranames) %>%
      pivot_longer(everything(),names_to = "para") %>%
      mutate(para=factor(para, levels=paranames))
  }
  
  df <- cbind(do.call("rbind", dfs),
              label=rep(labels, each=nsim*npara))
  
  p <- df %>%
    ggplot(aes(x=para, y=value, fill=label)) + 
    geom_boxplot(outlier.size=1) + 
    theme_bw(base_size = 14) +
    labs(y=ylab, title=title,fill="",x="") + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(),
          legend.margin = margin(t = -15, unit = "pt"),
          legend.position = "bottom") +
    scale_x_discrete(labels=parse(text=paranames[xbreaks]),
                     breaks=paranames[xbreaks])
  
  if (!is.null(ylim)) 
    p <- p + lims(y=ylim)
  
  p
}

# helper function
# using spec$d_pred in the global environment
get_df_v2 <- function(df){
  data.frame(cbind(spec$d_pred,df)) %>%
    `colnames<-`(c("x","y","value"))
}

# Obtain scale of the data to prepare for plots
# bl, bu, bby: for break points
# ll, uu: for limits
get_plot_scale_v3 <- function(rs,ni,j,dfs){
  
  values <- do.call("rbind",dfs)
  
  values <- values[,j]
  ll <- floor(min(values)*rs)/rs
  uu <- ceiling(max(values)*rs)/rs
  
  if (ll==uu){
    ll <- ll*9/10
    uu <- uu*11/10
  }
  
  bby <- signif((uu-ll)/ni, digits=1)
  
  list(ll=ll,uu=uu,bby=bby,bl=ll,bu=uu)
  
}

