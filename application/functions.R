library(matrixStats)

# adapted from v1, used for stan files where all hyperparameters are data input
runHSGP_v2 <- function(rs,ss,mod,i,seed=1,testshort=F,testfull=F){
  
  if (testshort){
    t_iter <- 200; t_burn <- 100
    refresh <- max(10, t_iter/10)
  } else {
    t_iter <- rs$niter; t_burn <- rs$nburn
    
    if (testfull){
      refresh <- max(10, t_iter/10)
    } else {
      refresh <- 0
    }
  }
  
  if (ss$isModels & !ss$isModeli){
    data <- list(n=ss$n,n_pred=ss$n_pred,p=ss$p,k=ss$k,m=ss$m,y=rs$y,s=rs$s,
                 delta=rs$delta,X=ss$X,W=ss$W,d=ss$d,d_pred=ss$d_pred,a=ss$a,
                 b=ss$b,kappa0=ss$kappa0,L=rs$c[[i]][1,,]*rs$S,M=rs$m[[i]][1,,],
                 a0=ss$a0,b0=ss$b0,a1=ss$a1,b1=ss$b1,al=ss$al,bl=ss$bl,tau=ss$tau)
    
    rs$fit[[i]] <- sampling(mod, data = data, init=0.5, chains = rs$nchain, 
                            warmup = t_burn,iter = t_iter, core=rs$nchain, 
                            thin = rs$nthin, seed = seed, include=T, refresh=refresh,
                            pars = c("l1","alpha1","beta","beta_w","s2","kappa",
                                     "lambda","theta_pred","alpha1_adj","loglik"))  
    
  } else if (ss$isModeli & ss$isModels){
    data <- list(n=ss$n,n_pred=ss$n_pred,p=ss$p,k=ss$k,y=rs$y,s=rs$s,
                 delta=rs$delta,X=ss$X,W=ss$W,d=ss$d,d_pred=ss$d_pred,a=ss$a,
                 b=ss$b,kappa0=ss$kappa0,L0=rs$c[[i]][1,,]*rs$S,M0=rs$m[[i]][1,,],
                 L1=rs$c[[i]][2,,]*rs$S,M1=rs$m[[i]][2,,],
                 a0=ss$a0,b0=ss$b0,a1=ss$a1,b1=ss$b1,al=ss$al,bl=ss$bl,tau=ss$tau)
    
    rs$fit[[i]] <- sampling(mod, data = data, init=0.5, chains = rs$nchain, 
                            warmup = t_burn,iter = t_iter, core=rs$nchain, 
                            thin = rs$nthin, seed = seed, include=T, refresh=refresh,
                            pars = c("l0","alpha0","l1","alpha1","beta","beta_w",
                                     "s2","kappa","lambda","theta_pred",
                                     "alpha0_adj","alpha1_adj","loglik"))      
    
  } else if (!ss$isModels & ss$isModeli){
    data <- list(n=ss$n,n_pred=ss$n_pred,p=ss$p,k=ss$k,y=rs$y,s=rs$s,
                 delta=rs$delta,X=ss$X,d=ss$d,d_pred=ss$d_pred,a=ss$a,
                 b=ss$b,kappa0=ss$kappa0,L0=rs$c[[i]][1,,]*rs$S,M0=rs$m[[i]][1,,],
                 a0=ss$a0,b0=ss$b0,a1=ss$a1,b1=ss$b1,al=ss$al,bl=ss$bl,tau=ss$tau)
    
    rs$fit[[i]] <- sampling(mod, data = data, init=0.5, chains = rs$nchain, 
                            warmup = t_burn,iter = t_iter, core=rs$nchain, 
                            thin = rs$nthin, seed = seed, include=T, refresh=refresh,
                            pars = c("l0","alpha0","beta","s2","kappa","lambda",
                                     "theta_pred","alpha0_adj","loglik")) 
  }
  
  return(rs)
}

updateHSGP_v1 <- function(rs,ss,i){
  S <- rs$S
  d <- 2
  
  if (ss$isModeli & ss$isModels) {
    nGP <- 2
    l_labels <- c("l0","l1")
  } else {
    nGP <- 1
    if (ss$isModeli) {
      l_labels <- "l0"
    } else if (ss$isModels) {
      l_labels <- "l1"
    }
  }
  
  m <- ss$m
  
  # we assume each GP has only one lengthscale parameter, for all the dimensions
  # but HSGP setup would require one lengthscale for each dimension during configuration
  # because the length of each dimension could be different
  
  ci <- mi <- lij <- array(0,dim=c(nGP,d,m))
  li <- li_hat <- matrix(0,nGP,m)
  diagi <- check <- matrix(F,nGP,m)
  
  if (i==1){
    li <- matrix(0.5*min(S),nGP,m)
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

# adapted from v10 above, try not to specify title="" if user input is NULL
get_plot1_v11 <- function(geo_plot,map_sf,color_plot,geo_obs=NULL,text_df=NULL,
                          title=NULL,legend=NULL,m_scale=NULL,center=0){
  
  # this function visualizes colors on a continuous scale
  # geo_plot: sf object, geometry (point coordinates) for the color plot
  # map_sf: sf object, geometry for the overall map
  # color_plot: colors for geo_plot
  # geo_obs: sf object, geometry (point coordinates) for observed data points, if null, no observation points will be plotted
  # text_df: text labels to be annotated
  # m_scale: user-specified scales, typically used if this plot needs to be on the same scale as another plot
  #         if not input, this function will get the scale for color_plot within the function
  # center: center of the scales
  
  ncolor <- 10
  if (!is.null(geo_obs)) n <- dim(geo_obs)[1]
  legend <- ifelse(is.null(legend), "", legend)
  
  if (is.null(m_scale)){
    # rounding scalar
    rs <- 10
    # number of intervals
    ni <- 4
    
    ll <- floor(min(color_plot)*rs)/rs
    ul <- ceiling(max(color_plot)*rs)/rs
    
    if (ll==ul){
      ll <- ll*9/10
      ul <- ul*11/10
    }
    
    lby <- signif((ul-ll)/ni, digits=1)
    
    m_scale <- list(ll=ll,ul=ul,lby=lby)
  }
  
  if (m_scale$ll<center & m_scale$ul>center){
    if (center-m_scale$ll < m_scale$ul-center) {
      ncolor1 <- max(round(((center-m_scale$ll)/(m_scale$ul-center))*ncolor,0),1)
      ncolor2 <- ncolor
    } else {
      ncolor1 <- ncolor
      ncolor2 <- max(round(((m_scale$ul-center)/(center-m_scale$ll))*ncolor,0),1)
    }
    
    color_scale <- c(tail(viridis(15, option = "G"),ncolor1), "white", 
                     rev(tail(viridis(15, option = "A"),ncolor2)))
    color_values <- scales::rescale(seq(m_scale$ll,m_scale$ul, length.out=length(color_scale)))
    
  } else {
    color_values <- c(0,1)
    if (m_scale$ll>=center){
      color_scale <- rev(tail(viridis(15, option = "A"),ncolor))
    } else if (m_scale$ul<=center){
      color_scale <- tail(viridis(15, option = "G"),ncolor)
    }
  }
  
  p <- ggplot()+
    geom_sf(data=geo_plot, aes(color=color_plot), shape=15, size=10) +
    geom_sf(data=map_sf, fill=NA, color="#bdbdbd", size=1) + 
    theme_bw() + 
    theme(panel.border = element_blank(),
          legend.title = element_text(size=8),
          legend.text = element_text(size=8),
          legend.position = "bottom",
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          legend.margin = margin(t = -15, unit = "pt")) + 
    scale_color_gradientn(colors=color_scale, values = color_values,
                          breaks=round(seq(m_scale$ll,m_scale$ul,by=m_scale$lby),2),
                          limits=c(m_scale$ll,m_scale$ul)) 
  
  if (is.null(title)){
    p <- p + labs(x="",y="",color=legend)
  } else {
    p <- p + labs(title=title,x="",y="",color=legend)
  }
  
  if (!is.null(geo_obs))
    p <- p + geom_sf(data=geo_obs, col = "#969696", shape=16, size=1)
  
  if (!is.null(text_df))
    p <- p + geom_text(data = text_df, aes(lon, lat, label = city), size = 3,
                       fontface="bold")
  
  p
}


# adapted from above, changed the calculation for theta1_pred_level to account for 
# differences between average then exp vs exp then average
extract_beta_v2 <- function(spec,fit){
  
  beta_names <- grep("beta", names(fit), value = TRUE)
  beta <- as.matrix(fit,pars=beta_names)
  n_psim <- dim(beta)[1]
  
  if (spec$isModels){
    nGP <- sum(spec$isModeli+spec$isModels)
    theta_pred <- array(rstan::extract(fit,"theta_pred")$theta_pred,
                        dim=c(n_psim,spec$m,spec$n_pred,nGP))
    idx <- which(str_detect(colnames(beta),"beta_w"))
    
    theta1_pred <- vapply(c(1:spec$n_pred), function(j) theta_pred[,,j,nGP]+beta[,idx], matrix(0,n_psim,m)) #npsim x m x n_pred
    
    theta1_pred_level <- matrix(0,n_psim,spec$m)
    for (j in 1:spec$m){
      theta1_pred_level[,j] <- log(rowMeans(exp(theta1_pred[,j,])))
    }
    beta[,idx] <- theta1_pred_level
  }
  
  return(beta)
}

# get the plotting scales for theta related figures
# this version takes in matrices instead of dataframes
get_plot_scale_v2 <- function(rs,ni,j,dfs){
  
  values <- do.call("rbind",dfs)
  
  values <- values[,j]
  ll <- floor(min(values)*rs)/rs
  ul <- ceiling(max(values)*rs)/rs
  
  if (ll==ul){
    ll <- ll*9/10
    ul <- ul*11/10
  }
  
  lby <- signif((ul-ll)/ni, digits=1)
  
  list(ll=ll,ul=ul,lby=lby)
  
}


# this version can work with different risk types
getHR_multi_v3 <- function(res,rowlabels,ordered_labels=NULL,p=NULL,risktype=1){
  
  idx <- which(str_detect(colnames(res),paste0("beta.*(,",risktype,"|\\[",risktype,"\\])")))
  if (is.null(p)) p <- length(idx)
  
  p_beta <- res[,idx]
  e_p_beta <- exp(p_beta)
  HR_m <- colMeans(e_p_beta)
  HR_u <- apply(e_p_beta,2,quantile,0.975)
  HR_l <-  apply(e_p_beta,2,quantile,0.025)
  
  df <- data.frame(HR_m=HR_m, HR_l=HR_l, HR_u=HR_u) %>%
    mutate(label=rowlabels) %>%
    arrange(HR_m)
  
  if (!is.null(ordered_labels)){
    df <- df[match(ordered_labels,df$label),]
  }
  
  return(df)
}


# adapted from plotHRs_v2, added option to change legend position
plotHRs_v3 <- function(dfs,model,ordered_labels=NULL,xlim=NULL,yaxt=T,cex=0.8,
                       legend_labels=NULL,legend_position="topleft"){
  
  nr <- length(dfs)
  p <- dim(dfs[[1]])[1]
  tol <- 0.00001
  
  # make sure all the tables have the same order
  if (is.null(ordered_labels))
    ordered_labels <- dfs[[1]]$label
  
  dfs <- lapply(dfs, function(x) x[match(ordered_labels,x$label),])
  
  if (is.null(xlim)){
    plot(range(dfs[[1]][,1:3]), c(1,p), type="n", yaxt="n", main=model,
         ylab="", xlab="Hazard Ratio",cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex)
  } else {
    plot(range(dfs[[1]][,1:3]), c(1,p), type="n", yaxt="n", main=model, xlim=xlim,
         ylab="", xlab="Hazard Ratio",cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex)
  }
  
  for (j in 1:nr){
    offset <- (j-1)*0.1
    for(i in 1:p){
      
      x <- dfs[[j]]$HR_m[i]
      
      if (!is.na(x)){
        
        if (dfs[[j]]$HR_l[i]>1 | dfs[[j]]$HR_u[i]<1){
          mycol <- "red"
        } else {
          mycol <- "black"
        }
        
        if (abs(dfs[[j]]$HR_l[i]-dfs[[j]]$HR_u[i])<tol)
          mycol <- "grey"
        
        points( dfs[[j]][i,2:3], rep(i,2)-offset, pch=124,cex=.6, col=mycol)
        points( x, i-offset, pch=16,cex=.6, col=mycol)
        segments( dfs[[j]][i,2], i-offset, dfs[[j]][i,3], i-offset, col=mycol, lty=j)        
      }
      
    }
  }
  
  abline(v=1, col="blue")
  if (yaxt)
    axis(2, at = c(1:p), labels=ordered_labels, las=2, cex.axis=cex)  
  
  if (!is.null(legend_labels))
    legend(legend_position, legend = legend_labels,lty = c(1:nrun),bty = "n") 
} 

get_df_run <- function(df,outcome){
  
  if (outcome == "HP"){
    df_run <- df %>%
      mutate(time = pmin(time_c,time_d,time_hp,na.rm=T),
             event = case_when(time==time_c ~ 0,
                               time==time_hp ~ 1,
                               time==time_d ~ 2))
  }else if (outcome == "FA"){
    df_run <- df %>%
      mutate(time = pmin(time_c,time_d,time_fa,na.rm=T),
             event = case_when(time==time_c ~ 0,
                               time==time_fa ~ 1,
                               time==time_d ~ 2))
    
  }else if (outcome == "FR"){
    df_run <- df %>%
      mutate(time = pmin(time_c,time_d,time_fr,na.rm=T),
             event = case_when(time==time_c ~ 0,
                               time==time_fr ~ 1,
                               time==time_d ~ 2))
    
  }else if (outcome == "ED"){
    df_run <- df %>%
      mutate(time = pmin(time_c,time_d,time_ed,na.rm=T),
             event = case_when(time==time_c ~ 0,
                               time==time_ed ~ 1,
                               time==time_d ~ 2))
    
  }
  
  df_run
  
}
