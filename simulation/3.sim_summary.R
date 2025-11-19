
specnum <- 1
runnum <- 1
nsim <- 500

# folder path which stores the "functions.R" file and the spec file, make sure the path ends with "/"
infolder <- ""
# folder path with GP run results, make sure the path ends with "/"
resfolder_GP <- ""
# folder path with HSGP run results, make sure the path ends with "/"
resfolder_HSGP <- ""
# folder path with frequentist methods run results, make sure the path ends with "/"
resfolder_freq <- ""

library(rstan)
library(abind)
library(cowplot)

source(paste0(infolder,"functions.R"))

# get summary statistics
res_GP <- get_stats_v5(infolder,resfolder_GP,infolder,"GP",specnum,runnum,nsim=nsim,n_cpus=6)
res_HSGP <- get_stats_v5(infolder, resfolder_HSGP,infolder,"HSGP",specnum,runnum,nsim=nsim,n_cpus=6)
res_freq <- get_stats_freq_v2(infolder, resfolder_freq,infolder,specnum,runnum,nsim=nsim)

# read in spec information and set up parameters
spec <- readRDS(paste0(infolder, "spec_",specnum,".rds"))
isModeli <- spec$isModeli
isModels <- spec$isModels
isJointwithin <- spec$isJointwithin
isModelrt <- spec$isModelrt

n <- spec$n; n_pred <- spec$n_pred; p <- spec$p; m <- spec$m

if (spec$pl){
  s <- spec$s
  tmp <- map(res_HSGP,`[[`,"lambda_std_diff")
  k <- min(sapply(tmp, function(x) dim(x)[1]))
} else {
  k <- spec$k
}

para_res <- get_labels_v2(spec)
paras <- para_res$para_latex
npara <- length(paras)

lambda_label <- matrix(nrow=k,ncol=m)
for (j in 1:m) {
  lambda_label[, j] <- sprintf('lambda["%d,%d"]', j, 1:k)
}
lambdabreaks <- c(1,seq(5,k,by=5))
GP_label <- c("BSp GP","BSp HSGP")

# rounding scalar
rs <- 10
# number of intervals
ni <- 4

df_d <- data.frame(spec$d) %>%
  `colnames<-`(c("x","y"))

if (spec$isModeli)
  f0_true <- t(spec$f_true_check[,1,])

if (spec$isModels)
  f1_true <- t(spec$f_true_check[,spec$s_idx,])

# compute summary statistics
if (isModeli){
  GP_theta0_m <- apply(abind( map(res_GP,`[[`,"theta0_m"), along=3 ),2,rowMeans)
  GP_theta0_std_diff <- apply(abind( map(res_GP,`[[`,"theta0_std_diff"), along=3 ),2,rowMeans)
  GP_theta0_coverage <- apply(abind( map(res_GP,`[[`,"theta0_coverage"), along=3 ),2,rowMeans)
  GP_theta0_sd <- apply(abind( map(res_GP,`[[`,"theta0_sd"), along=3 ),2,rowMeans)
  GP_theta0_rmse <- apply(abind( map(res_GP,`[[`,"theta0_rmse"), along=3 ),2,rowMeans)
  
  HSGP_theta0_m <- apply(abind( map(res_HSGP,`[[`,"theta0_m"), along=3 ),2,rowMeans)
  HSGP_theta0_std_diff <- apply(abind( map(res_HSGP,`[[`,"theta0_std_diff"), along=3 ),2,rowMeans)
  HSGP_theta0_coverage <- apply(abind( map(res_HSGP,`[[`,"theta0_coverage"), along=3 ),2,rowMeans)
  HSGP_theta0_sd <- apply(abind( map(res_HSGP,`[[`,"theta0_sd"), along=3 ),2,rowMeans)
  HSGP_theta0_rmse <- apply(abind( map(res_HSGP,`[[`,"theta0_rmse"), along=3 ),2,rowMeans)
}

if (isModels){
  GP_theta1_m <- apply(abind( map(res_GP,`[[`,"theta1_m"), along=3 ),2,rowMeans)
  GP_theta1_std_diff <- apply(abind( map(res_GP,`[[`,"theta1_std_diff"), along=3 ),2,rowMeans)
  GP_theta1_coverage <- apply(abind( map(res_GP,`[[`,"theta1_coverage"), along=3 ),2,rowMeans)
  GP_theta1_sd <- apply(abind( map(res_GP,`[[`,"theta1_sd"), along=3 ),2,rowMeans)
  GP_theta1_rmse <- apply(abind( map(res_GP,`[[`,"theta1_rmse"), along=3 ),2,rowMeans)
  
  HSGP_theta1_m <- apply(abind( map(res_HSGP,`[[`,"theta1_m"), along=3 ),2,rowMeans)
  HSGP_theta1_std_diff <- apply(abind( map(res_HSGP,`[[`,"theta1_std_diff"), along=3 ),2,rowMeans)
  HSGP_theta1_coverage <- apply(abind( map(res_HSGP,`[[`,"theta1_coverage"), along=3 ),2,rowMeans)
  HSGP_theta1_sd <- apply(abind( map(res_HSGP,`[[`,"theta1_sd"), along=3 ),2,rowMeans)
  HSGP_theta1_rmse <- apply(abind( map(res_HSGP,`[[`,"theta1_rmse"), along=3 ),2,rowMeans)
  
}

GP_ESS <- do.call("rbind",map(res_GP,`[[`,"ESS"))
GP_RT <- do.call("rbind",map(res_GP,`[[`,"RT"))
GP_para_std_diff <- do.call("rbind",map(res_GP,`[[`,"para_std_diff"))
GP_para_m_diff <- do.call("rbind",map(res_GP,`[[`,"para_m_diff"))
GP_para_m_rmse <- sqrt(colMeans(GP_para_m_diff^2))
GP_para_rmse <- do.call("rbind",map(res_GP,`[[`,"para_rmse"))
GP_para_coverage <- colMeans(do.call("rbind",map(res_GP,`[[`,"para_coverage")))
if (spec$pl){
  GP_lambda_m_list <- lapply(map(res_GP,`[[`,"lambda_m"), function(x) x[1:k,])
  GP_lambda_std_diff_list <- lapply(map(res_GP,`[[`,"lambda_std_diff"), function(x) x[1:k,])
  GP_lambda_rmse_list <- lapply(map(res_GP,`[[`,"lambda_rmse"), function(x) x[1:k,])
  GP_lambda_coverage_list <- lapply(map(res_GP,`[[`,"lambda_coverage"), function(x) x[1:k,])
  GP_lambda_m <- abind(GP_lambda_m_list, along=3)
  GP_lambda_std_diff <- abind(GP_lambda_std_diff_list, along=3)
  GP_lambda_rmse <- abind(GP_lambda_rmse_list, along=3)
  GP_lambda_coverage <- apply(abind( GP_lambda_coverage_list, along=3 ), c(1,2), mean)
} else {
  sim_s <- do.call(rbind,map(res_GP,`[[`,"s"))
  GP_lambda_m <- abind( map(res_GP,`[[`,"lambda_m"), along=3 )
  GP_lambda_std_diff <- abind( map(res_GP,`[[`,"lambda_std_diff"), along=3 )
  GP_lambda_rmse <- abind( map(res_GP,`[[`,"lambda_rmse"), along=3 )
  GP_lambda_coverage <- apply(abind( map(res_GP,`[[`,"lambda_coverage"), along=3 ), c(1,2), mean)
}

alpha_idx <- str_detect(paras,"alpha")
GP_alpha_std_diff <- GP_para_std_diff[,alpha_idx]
GP_alpha_rmse <- GP_para_rmse[,alpha_idx]
GP_alpha_coverage <- GP_para_coverage[alpha_idx]

HSGP_ESS <- do.call("rbind",map(res_HSGP,`[[`,"ESS"))
HSGP_RT <- do.call("rbind",map(res_HSGP,`[[`,"RT"))
HSGP_para_std_diff <- do.call("rbind",map(res_HSGP,`[[`,"para_std_diff"))
HSGP_para_m_diff <- do.call("rbind",map(res_HSGP,`[[`,"para_m_diff"))
HSGP_para_m_rmse <- sqrt(colMeans(HSGP_para_m_diff^2))
HSGP_para_rmse <- do.call("rbind",map(res_HSGP,`[[`,"para_rmse"))
HSGP_para_coverage <- colMeans(do.call("rbind",map(res_HSGP,`[[`,"para_coverage")))
HSGP_alpha_adj_std_diff <- do.call("rbind",map(res_HSGP,`[[`,"alpha_adj_std_diff"))
HSGP_alpha_adj_rmse <- do.call("rbind",map(res_HSGP,`[[`,"alpha_adj_rmse"))
HSGP_alpha_adj_coverage <- colMeans(do.call("rbind",map(res_HSGP,`[[`,"alpha_adj_coverage")))
if (spec$pl){
  HSGP_lambda_m_list <- lapply(map(res_HSGP,`[[`,"lambda_m"), function(x) x[1:k,])
  HSGP_lambda_std_diff_list <- lapply(map(res_HSGP,`[[`,"lambda_std_diff"), function(x) x[1:k,])
  HSGP_lambda_rmse_list <- lapply(map(res_HSGP,`[[`,"lambda_rmse"), function(x) x[1:k,])
  HSGP_lambda_coverage_list <- lapply(map(res_HSGP,`[[`,"lambda_coverage"), function(x) x[1:k,])
  HSGP_lambda_m <- abind(HSGP_lambda_m_list, along=3)
  HSGP_lambda_std_diff <- abind(HSGP_lambda_std_diff_list, along=3)
  HSGP_lambda_rmse <- abind(HSGP_lambda_rmse_list, along=3)
  HSGP_lambda_coverage <- apply(abind( HSGP_lambda_coverage_list, along=3 ), c(1,2), mean)
} else {
  # sim_s <- do.call(rbind,map(res_HSGP,`[[`,"s"))
  HSGP_lambda_m <- abind( map(res_HSGP,`[[`,"lambda_m"), along=3 )
  HSGP_lambda_std_diff <- abind( map(res_HSGP,`[[`,"lambda_std_diff"), along=3 )
  HSGP_lambda_rmse <- abind( map(res_HSGP,`[[`,"lambda_rmse"), along=3 )
  HSGP_lambda_coverage <- apply(abind( map(res_HSGP,`[[`,"lambda_coverage"), along=3 ), c(1,2), mean)
}

coxph1_para_std_diff <- do.call("rbind",map(res_freq,`[[`,"coxph1_std_diff"))
coxph2_para_std_diff <- do.call("rbind",map(res_freq,`[[`,"coxph2_std_diff"))
coxph1_para_m_diff <- do.call("rbind",map(res_freq,`[[`,"coxph1_m_diff"))
coxph2_para_m_diff <- do.call("rbind",map(res_freq,`[[`,"coxph2_m_diff"))
coxph1_rmse <- sqrt(colMeans(coxph1_para_m_diff^2, na.rm=T))
coxph2_rmse <- sqrt(colMeans(coxph2_para_m_diff^2,na.rm=T))

### Plotting

# Figure 5 and Figure 6 in the article

# bl, bu, bby: lower bound, upper bound and "by" for breaks
# ll, uu: lower and upper bounds for limits
myscale <- list(bl=-1.2,bu=1.2,bby=0.6,ll=-1.3,uu=1.3)
pA <- to_get_plot_v5(get_plot1_v6,"A",rt=1,myscale=myscale,
                     c("Truth",GP_label), list(f0_true,GP_theta0_m,HSGP_theta0_m),
                     legend="Spaital intercept")

myscale <- list(bl=-0.6,bu=1.2,bby=0.6,ll=-0.6,uu=1.3)
pB <- to_get_plot_v5(get_plot1_v6,"B",rt=1,myscale=myscale,
                     c("Truth",GP_label), list(f1_true,GP_theta1_m,HSGP_theta1_m),
                     legend="Spatial slope")

pC <- to_get_plot_v5(get_plot1_v6,"A",rt=1, GP_label,
                     list(GP_theta0_sd,HSGP_theta0_sd), legend="Spatial intercept posterior SD")

myscale <- list(bl=0.4,bu=0.7,bby=0.1,ll=0.4,uu=0.7)
pD <- to_get_plot_v5(get_plot1_v6,"B",rt=1, GP_label,myscale=myscale,
                     list(GP_theta1_sd,HSGP_theta1_sd), legend="Spatial slope posterior SD")


F5 <- plot_grid(
  pA,
  NULL,
  pB,
  ncol = 1,
  rel_heights = c(1, 0.05,1)
)

F6 <- plot_grid(
  pC, pD,
  nrow = 1,
  rel_widths = c(1, 1)
)

# Figure 3 in the supplement
myscale <- list(bl=-1,bu=1,bby=0.5,ll=-1.1,uu=1)
pA <- to_get_plot_v5(get_plot1_v6,"A",rt=2,myscale=myscale,
                     c("Truth",GP_label), list(f0_true,GP_theta0_m,HSGP_theta0_m),
                     legend="Spaital intercept")

myscale <- list(bl=-2,bu=0.5,bby=0.5,ll=-2,uu=0.8)
pB <- to_get_plot_v5(get_plot1_v6,"B",rt=2,myscale=myscale,
                     c("Truth",GP_label), list(f1_true,GP_theta1_m,HSGP_theta1_m),
                     legend="Spatial slope")

S <- plot_grid(
  pA,
  NULL,
  pB,
  ncol = 1,
  rel_heights = c(1, 0.05,1)
)

# Figure 4 in the supplement
pA <- to_get_plot_v5(get_plot1_v6,"A",rt=2, GP_label,
                     list(GP_theta0_sd,HSGP_theta0_sd), legend="Spatial intercept posterior SD")

pB <- to_get_plot_v5(get_plot1_v6,"B",rt=2, GP_label,
                     list(GP_theta1_sd,HSGP_theta1_sd), legend="Spatial slope posterior SD")

S <- plot_grid(
  pA, pB,
  nrow = 1,
  rel_widths = c(1, 1)
)

# Figure 5 in the supplement

pA <- to_get_plot_v5(get_plot1_v6,"A",rt=1,GP_label,
                     list(GP_theta0_rmse,HSGP_theta0_rmse), legend="Spatial intercept RMSE")
pB <- to_get_plot_v5(get_plot1_v6,"B",rt=1,GP_label,
                     list(GP_theta1_rmse,HSGP_theta1_rmse), legend="Spatial slope RMSE")

pC <- to_get_plot_v5(get_plot1_v6,"C",rt=2,GP_label,
                     list(GP_theta0_rmse,HSGP_theta0_rmse), legend="Spatial intercept RMSE")
pD <- to_get_plot_v5(get_plot1_v6,"D",rt=2,GP_label,
                     list(GP_theta1_rmse,HSGP_theta1_rmse), legend="Spatial slope RMSE")

S <- plot_grid(
  pA,pB,pC,pD,
  nrow=2, ncol=2,
  rel_widths = c(1,1),
  rel_heights = c(1,1)
)

# Figure 6 in the supplement
pA <- to_get_plot_v5(get_plot2_v3,"A",rt=1,GP_label,
                     list(GP_theta0_coverage,HSGP_theta0_coverage), legend="Spatial intercept posterior coverage")
pB <- to_get_plot_v5(get_plot2_v3,"B",rt=1,GP_label,
                     list(GP_theta1_coverage,HSGP_theta1_coverage), legend="Spatial slope posterior coverage")
pC <- to_get_plot_v5(get_plot2_v3,"C",rt=2,GP_label,
                     list(GP_theta0_coverage,HSGP_theta0_coverage), legend="Spatial intercept posterior coverage")
pD <- to_get_plot_v5(get_plot2_v3,"D",rt=2,GP_label,
                     list(GP_theta1_coverage,HSGP_theta1_coverage), legend="Spatial slope posterior coverage")

S <- plot_grid(
  pA,pB,pC,pD,
  nrow=2, ncol=2,
  rel_widths = c(1,1),
  rel_heights = c(1,1)
)

pdf(file = paste0(outfolder,"S_sim_sp_cov_v3.pdf"), width=10, height=6)
S
dev.off()

# Figure 4 in the article

lambda_f <- partial(spec$lambda_f_check,spec=spec)
vec_lambda_f <- Vectorize(lambda_f)

min_s <- min(sim_s[,k+1])

j <- 1

# deduct 0.1 from min_s to avoid out of bounds error below while using findInterval
p_GP <- p_HSGP <- data.frame(x = 0:(min_s-0.1)) %>% 
  ggplot(aes(x = x))

for (i in 1:nsim){
  
  GP_lambda_ij <- function(t,i,j){
    idx <- findInterval(t,sim_s[i,])
    GP_lambda_m[idx,j,i]
  }
  HSGP_lambda_ij <- function(t,i,j){
    idx <- findInterval(t,sim_s[i,])
    HSGP_lambda_m[idx,j,i]        
  }
  
  p_GP <- p_GP + stat_function(fun = GP_lambda_ij, args=list(i=i,j=j), col="red", 
                               alpha=0.5)
  p_HSGP <- p_HSGP + stat_function(fun = HSGP_lambda_ij, args=list(i=i,j=j), col="red", 
                                   alpha=0.5)
  
}

p_GP <- p_GP +
  stat_function(fun = vec_lambda_f,args=list(j=j)) +
  theme_bw() +
  labs(x="Time since baseline",y="Baseline hazard rate",title="BSp GP") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p_HSGP <- p_HSGP +
  stat_function(fun = vec_lambda_f,args=list(j=j)) +
  theme_bw() +
  labs(x="Time since baseline",y="Baseline hazard rate",title="BSp HSGP") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())

F4 <- plot_grid(
  p_GP, p_HSGP,
  labels = c("A", "B"),
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1)
)


# Figure 2 in the supplement
j <- 2

# deduct 0.1 from min_s to avoid out of bounds error below while using findInterval
p_GP <- p_HSGP <- data.frame(x = 0:(min_s-0.1)) %>% 
  ggplot(aes(x = x))

for (i in 1:nsim){
  
  GP_lambda_ij <- function(t,i,j){
    idx <- findInterval(t,sim_s[i,])
    GP_lambda_m[idx,j,i]
  }
  HSGP_lambda_ij <- function(t,i,j){
    idx <- findInterval(t,sim_s[i,])
    HSGP_lambda_m[idx,j,i]        
  }
  
  p_GP <- p_GP + stat_function(fun = GP_lambda_ij, args=list(i=i,j=j), col="red", 
                               alpha=0.5)
  p_HSGP <- p_HSGP + stat_function(fun = HSGP_lambda_ij, args=list(i=i,j=j), col="red", 
                                   alpha=0.5)
  
}

p_GP <- p_GP +
  stat_function(fun = vec_lambda_f,args=list(j=j)) +
  theme_bw() +
  labs(x="Time since baseline",y="Baseline hazard rate",title="BSp GP") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p_HSGP <- p_HSGP +
  stat_function(fun = vec_lambda_f,args=list(j=j)) +
  theme_bw() +
  labs(x="Time since baseline",y="Baseline hazard rate",title="BSp HSGP") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())

F4 <- plot_grid(
  p_GP, p_HSGP,
  labels = c("A", "B"),
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1)
)

# Figure 1 in the supplement
ylim <- c(0,2)

j <- 1
pA <- to_get_plot3_v2(list(t(GP_lambda_rmse[,j,]),t(HSGP_lambda_rmse[,j,])),
                      lambda_label[,j],GP_label,paste0("Risk type ",j),
                      ylim,xbreaks=lambdabreaks,ylab="RMSE")
j <- 2 
pB <- to_get_plot3_v2(list(t(GP_lambda_rmse[,j,]),t(HSGP_lambda_rmse[,j,])),
                      lambda_label[,j],GP_label,paste0("Risk type ",j),
                      ylim,xbreaks=lambdabreaks,ylab="RMSE")

S <- plot_grid(
  pA, pB,
  labels = c("A", "B"),
  label_fontface = "bold",
  nrow = 2
)


