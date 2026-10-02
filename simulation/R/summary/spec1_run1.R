specnum <- 1
runnum <- 1
FQ_runnum <- 1

# 20260712 notes
# 1. simplify2array is faster than abind

nsim <- 1000

folder <- "simulation/"
source(paste0(folder, "R/functions.R"))
outfolder <- paste0(folder,"fig/")
if (!dir.exists(outfolder)) dir.create(outfolder, recursive = TRUE)

library(ggpubr)
library(cowplot)
library(knitr)

get_stats_v8(folder,"HSGP",specnum,runnum,nsim=1000,n_cpus=40)
get_stats_v8(folder,"GP",specnum,runnum,nsim=1000,n_cpus=40)
get_stats_freq_v3(folder,specnum,FQ_runnum,nsim=1000,n_cpus=8)

res_HSGP <- readRDS(paste0(folder,"res/stats/HSGP_spec",specnum,"_run",runnum,"_stats_v8.rds"))
res_GP <- readRDS(paste0(folder,"res/stats/GP_spec",specnum,"_run",runnum,"_stats_v8.rds"))
res_FQ <- readRDS(paste0(folder,"res/stats/FQ_spec",specnum,"_run",FQ_runnum,"_stats_v3.rds"))

spec <- readRDS(paste0(folder,"spec/spec", specnum,".rds"))
isModeli <- spec$isModeli
isModels <- spec$isModels
n <- spec$n; n_pred <- spec$n_pred; p <- spec$p; m <- spec$m

if (spec$pl){
  s <- spec$s
  tmp <- map(res_HSGP,`[[`,"lambda_m_diff")
  k <- min(sapply(tmp, function(x) dim(x)[1]))
} else {
  k <- spec$k
}

para_res <- get_labels_v3(spec,"GP")
paras <- para_res$para_latex
npara <- length(paras)

df_d <- data.frame(spec$d) %>%
  `colnames<-`(c("x","y"))

t_cif <- spec$t_cif; xs <- spec$cif_xs
nxs <- length(xs)

# choose to present year 5 CIF
cif_i <- 5
  
if (spec$isModeli)
  f0_true <- t(spec$f_true_check[,1,])

if (spec$isModels)
  f1_true <- t(spec$f_true_check[,spec$s_idx,])

# rounding scalar
rs <- 10
# number of intervals
ni <- 4

GP_label <- c("BSp GP","BSp HSGP")

# check that GP and HSGP runs are using the same data
# and that HSGP runs have converged
# check_y(folder,"GP","HSGP",specnum,runnum,runnum,nsim,n_cpus=2)


# n_pred x nxs x m (checked correct)
GP_cif_m <- apply(simplify2array(map(res_GP, `[[`, "cif_m")), MARGIN = c(2,3), rowMeans)
GP_cif_m_diff <- simplify2array(map(res_GP, `[[`, "cif_m_diff"))
GP_cif_coverage <- apply(simplify2array(map(res_GP, `[[`, "cif_coverage")), MARGIN = c(2,3), rowMeans)
GP_cif_sd <- apply(simplify2array(map(res_GP, `[[`, "cif_sd")), MARGIN = c(2,3), rowMeans)
GP_cif_rmse <- apply(simplify2array(map(res_GP, `[[`, "cif_rmse")), MARGIN = c(2,3), rowMeans)
GP_cif_m_rmse <- sqrt(apply(GP_cif_m_diff^2, MARGIN = c(1,2,3), FUN = mean))

HSGP_cif_m <- apply(simplify2array(map(res_HSGP, `[[`, "cif_m")), MARGIN = c(2,3), rowMeans)
HSGP_cif_m_diff <- simplify2array(map(res_HSGP, `[[`, "cif_m_diff"))
HSGP_cif_coverage <- apply(simplify2array(map(res_HSGP, `[[`, "cif_coverage")), MARGIN = c(2,3), rowMeans)
HSGP_cif_sd <- apply(simplify2array(map(res_HSGP, `[[`, "cif_sd")), MARGIN = c(2,3), rowMeans)
HSGP_cif_rmse <- apply(simplify2array(map(res_HSGP, `[[`, "cif_rmse")), MARGIN = c(2,3), rowMeans)
HSGP_cif_m_rmse <- sqrt(apply(HSGP_cif_m_diff^2, MARGIN = c(1,2,3), FUN = mean))

if (isModeli){

  GP_theta0_m <- apply(simplify2array( map(res_GP,`[[`,"theta0_m")),2,rowMeans)
  GP_theta0_m_diff <- simplify2array(map(res_GP, `[[`, "theta0_m_diff"))
  GP_theta0_coverage <- apply(simplify2array( map(res_GP,`[[`,"theta0_coverage")),2,rowMeans)
  GP_theta0_sd <- apply(simplify2array( map(res_GP,`[[`,"theta0_sd")),2,rowMeans)
  GP_theta0_rmse <- apply(simplify2array( map(res_GP,`[[`,"theta0_rmse")),2,rowMeans)
  GP_theta0_m_rmse <- sqrt(apply(GP_theta0_m_diff^2, MARGIN = c(1,2), FUN = mean))
  
  HSGP_theta0_m <- apply(simplify2array( map(res_HSGP,`[[`,"theta0_m")),2,rowMeans)
  HSGP_theta0_m_diff <- simplify2array(map(res_HSGP, `[[`, "theta0_m_diff"))
  HSGP_theta0_coverage <- apply(simplify2array( map(res_HSGP,`[[`,"theta0_coverage")),2,rowMeans)
  HSGP_theta0_sd <- apply(simplify2array( map(res_HSGP,`[[`,"theta0_sd")),2,rowMeans)
  HSGP_theta0_rmse <- apply(simplify2array( map(res_HSGP,`[[`,"theta0_rmse")),2,rowMeans)
  HSGP_theta0_m_rmse <- sqrt(apply(HSGP_theta0_m_diff^2, MARGIN = c(1,2), FUN = mean))
}

if (isModels){
  GP_theta1_m <- apply(simplify2array( map(res_GP,`[[`,"theta1_m")),2,rowMeans)
  GP_theta1_m_diff <- simplify2array(map(res_GP, `[[`, "theta1_m_diff"))
  GP_theta1_coverage <- apply(simplify2array( map(res_GP,`[[`,"theta1_coverage")),2,rowMeans)
  GP_theta1_sd <- apply(simplify2array( map(res_GP,`[[`,"theta1_sd")),2,rowMeans)
  GP_theta1_rmse <- apply(simplify2array( map(res_GP,`[[`,"theta1_rmse")),2,rowMeans)
  GP_theta1_m_rmse <- sqrt(apply(GP_theta1_m_diff^2, MARGIN = c(1,2), FUN = mean))
  
  HSGP_theta1_m <- apply(simplify2array( map(res_HSGP,`[[`,"theta1_m")),2,rowMeans)
  HSGP_theta1_m_diff <- simplify2array(map(res_HSGP, `[[`, "theta1_m_diff"))
  HSGP_theta1_coverage <- apply(simplify2array( map(res_HSGP,`[[`,"theta1_coverage")),2,rowMeans)
  HSGP_theta1_sd <- apply(simplify2array( map(res_HSGP,`[[`,"theta1_sd")),2,rowMeans)
  HSGP_theta1_rmse <- apply(simplify2array( map(res_HSGP,`[[`,"theta1_rmse")),2,rowMeans)
  HSGP_theta1_m_rmse <- sqrt(apply(HSGP_theta1_m_diff^2, MARGIN = c(1,2), FUN = mean))
  
}

GP_ESS <- do.call("rbind",map(res_GP,`[[`,"ESS"))
GP_RT <- do.call("rbind",map(res_GP,`[[`,"RT"))
GP_para_m_diff <- do.call("rbind",map(res_GP,`[[`,"para_m_diff"))
GP_para_m_rmse <- sqrt(colMeans(GP_para_m_diff^2))
GP_para_rmse <- do.call("rbind",map(res_GP,`[[`,"para_rmse"))
GP_para_coverage <- colMeans(do.call("rbind",map(res_GP,`[[`,"para_coverage")))
if (spec$pl){
  GP_lambda_m_list <- lapply(map(res_GP,`[[`,"lambda_m"), function(x) x[1:k,])
  GP_lambda_m_diff_list <- lapply(map(res_GP,`[[`,"lambda_m_diff"), function(x) x[1:k,])
  GP_lambda_rmse_list <- lapply(map(res_GP,`[[`,"lambda_rmse"), function(x) x[1:k,])
  GP_lambda_coverage_list <- lapply(map(res_GP,`[[`,"lambda_coverage"), function(x) x[1:k,])
  GP_lambda_m <- simplify2array(GP_lambda_m_list)
  GP_lambda_m_diff <- simplify2array(GP_lambda_m_diff_list)
  GP_lambda_rmse <- simplify2array(GP_lambda_rmse_list)
  GP_lambda_m_rmse <- sqrt(apply(GP_lambda_m_diff^2, MARGIN = c(1,2), FUN = mean))
  GP_lambda_coverage <- apply(simplify2array( GP_lambda_coverage_list ), c(1,2), mean)
} else {
  sim_s <- do.call(rbind,map(res_GP,`[[`,"s"))
  GP_lambda_m <- simplify2array( map(res_GP,`[[`,"lambda_m") )
  GP_lambda_m_diff <- simplify2array( map(res_GP,`[[`,"lambda_m_diff") )
  GP_lambda_rmse <- simplify2array( map(res_GP,`[[`,"lambda_rmse") )
  GP_lambda_m_rmse <- sqrt(apply(GP_lambda_m_diff^2, MARGIN = c(1,2), FUN = mean))
  GP_lambda_coverage <- apply(simplify2array( map(res_GP,`[[`,"lambda_coverage") ), c(1,2), mean)
}

HSGP_ESS <- do.call("rbind",map(res_HSGP,`[[`,"ESS"))
HSGP_RT <- do.call("rbind",map(res_HSGP,`[[`,"RT"))
HSGP_para_m_diff <- do.call("rbind",map(res_HSGP,`[[`,"para_m_diff"))
HSGP_para_m_rmse <- sqrt(colMeans(HSGP_para_m_diff^2))
HSGP_para_rmse <- do.call("rbind",map(res_HSGP,`[[`,"para_rmse"))
HSGP_para_coverage <- colMeans(do.call("rbind",map(res_HSGP,`[[`,"para_coverage")))
if (spec$pl){
  HSGP_lambda_m_list <- lapply(map(res_HSGP,`[[`,"lambda_m"), function(x) x[1:k,])
  HSGP_lambda_m_diff_list <- lapply(map(res_HSGP,`[[`,"lambda_m_diff"), function(x) x[1:k,])
  HSGP_lambda_rmse_list <- lapply(map(res_HSGP,`[[`,"lambda_rmse"), function(x) x[1:k,])
  HSGP_lambda_coverage_list <- lapply(map(res_HSGP,`[[`,"lambda_coverage"), function(x) x[1:k,])
  HSGP_lambda_m <- simplify2array(HSGP_lambda_m_list)
  HSGP_lambda_m_diff <- simplify2array(HSGP_lambda_m_diff_list)
  HSGP_lambda_rmse <- simplify2array(HSGP_lambda_rmse_list)
  HSGP_lambda_m_rmse <- sqrt(apply(HSGP_lambda_m_diff^2, MARGIN = c(1,2), FUN = mean))
  HSGP_lambda_coverage <- apply(simplify2array( HSGP_lambda_coverage_list ), c(1,2), mean)
} else {
  # sim_s <- do.call(rbind,map(res_HSGP,`[[`,"s"))
  HSGP_lambda_m <- simplify2array( map(res_HSGP,`[[`,"lambda_m") )
  HSGP_lambda_m_diff <- simplify2array( map(res_HSGP,`[[`,"lambda_m_diff") )
  HSGP_lambda_rmse <- simplify2array( map(res_HSGP,`[[`,"lambda_rmse") )
  HSGP_lambda_m_rmse <- sqrt(apply(HSGP_lambda_m_diff^2, MARGIN = c(1,2), FUN = mean))
  HSGP_lambda_coverage <- apply(simplify2array( map(res_HSGP,`[[`,"lambda_coverage") ), c(1,2), mean)
}

coxph1_para_std_diff <- do.call("rbind",map(res_FQ,`[[`,"coxph1_std_diff"))
coxph2_para_std_diff <- do.call("rbind",map(res_FQ,`[[`,"coxph2_std_diff"))
coxph1_para_m_diff <- do.call("rbind",map(res_FQ,`[[`,"coxph1_m_diff"))
coxph2_para_m_diff <- do.call("rbind",map(res_FQ,`[[`,"coxph2_m_diff"))
coxph1_m_rmse <- sqrt(colMeans(coxph1_para_m_diff^2, na.rm=T))
coxph2_m_rmse <- sqrt(colMeans(coxph2_para_m_diff^2, na.rm=T))

# CIF and spatial slopes plot
myscale <- list(bl=0,bu=0.09,bby=0.03,ll=0,uu=0.09)
pA <- to_get_plot_v5(get_plot1_v6,"A",rt=1,myscale=myscale,
                     c("Truth",GP_label), list(t_cif[,cif_i,],GP_cif_m[,cif_i,],HSGP_cif_m[,cif_i,]),
                     legend="CIF at Year 5")

myscale <- list(bl=-0.8,bu=1.6,bby=0.8,ll=-0.8,uu=2.2)
pB <- to_get_plot_v5(get_plot1_v6,"B",rt=1,myscale=myscale,
                     c("Truth",GP_label), list(f1_true,GP_theta1_m,HSGP_theta1_m),
                     legend="Spatial slope")

# tried setting upper bound to 0.07, but the pattern is similar
myscale <- list(bl=0,bu=0.09,bby=0.03,ll=0,uu=0.09)
pC <- to_get_plot_v5(get_plot1_v6,"A",rt=1, GP_label,myscale=myscale,
                     list(GP_cif_sd[,cif_i,],HSGP_cif_sd[,cif_i,]), 
                     legend="CIF at Year 5 SD")

pD <- to_get_plot_v5(get_plot1_v6,"B",rt=1, GP_label,
                     list(GP_theta1_sd,HSGP_theta1_sd), legend="Spatial slope posterior SD")

row3 <- plot_grid(
  pC, pD,
  nrow = 1,
  rel_widths = c(1, 1)
)

F4 <- plot_grid(
  pA,
  NULL,
  pB,
  ncol = 1,
  rel_heights = c(1, 0.05,1)
)

pdf(file = paste0(outfolder,"F4_sim_sp_m_r1_v2.pdf"), width=10, height=8)
F4
dev.off()

pdf(file = paste0(outfolder,"F5_sim_sp_sd_r1_v2.pdf"), width=10, height=3)
row3
dev.off()

# Risk type 2 CIF and spatial slopes posterior mean
pA <- to_get_plot_v5(get_plot1_v6,"A",rt=2,
                     c("Truth",GP_label), list(t_cif[,cif_i,],GP_cif_m[,cif_i,],HSGP_cif_m[,cif_i,]),
                     legend="CIF at Year 5")

myscale <- list(bl=-3,bu=1,bby=1,ll=-3,uu=1)
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

pdf(file = paste0(outfolder,"S_sim_sp_m_r2_v2.pdf"), width=10, height=8)
S
dev.off()

# Risk type 2 CIF and posterior slope posterior SD
pA <- to_get_plot_v5(get_plot1_v6,"A",rt=2, GP_label,
                     list(GP_cif_sd[,cif_i,],HSGP_cif_sd[,cif_i,]), legend="CIF at Year 5 SD")

pB <- to_get_plot_v5(get_plot1_v6,"B",rt=2, GP_label,
                     list(GP_theta1_sd,HSGP_theta1_sd), legend="Spatial slope posterior SD")

S <- plot_grid(
  pA, pB,
  nrow = 1,
  rel_widths = c(1, 1)
)

pdf(file = paste0(outfolder,"S_sim_sp_sd_r2_v2.pdf"), width=10, height=3)
S
dev.off()

# Spatial intercepts

myscale <- list(bl=-2,bu=1,bby=1,ll=-2.2,uu=1.3)
pA <- to_get_plot_v5(get_plot1_v6,"A",rt=1,myscale=myscale,
                     c("Truth",GP_label), list(f0_true,GP_theta0_m,HSGP_theta0_m),
                     legend="Spaital intercept")
myscale <- list(bl=-2,bu=1,bby=1,ll=-2,uu=1.1)
pB <- to_get_plot_v5(get_plot1_v6,"B",rt=2,
                     c("Truth",GP_label), list(f0_true,GP_theta0_m,HSGP_theta0_m),
                     legend="Spaital intercept")
pC <- to_get_plot_v5(get_plot1_v6,"A",rt=1, GP_label,
                     list(GP_theta0_sd,HSGP_theta0_sd), legend="Spatial intercept posterior SD")
pD <- to_get_plot_v5(get_plot1_v6,"B",rt=2, GP_label,
                     list(GP_theta0_sd,HSGP_theta0_sd), legend="Spatial intercept posterior SD")

row3 <- plot_grid(
  pC, pD,
  nrow = 1,
  rel_widths = c(1, 1)
)

S <- plot_grid(
  pA,
  NULL,
  pB,
  ncol = 1,
  rel_heights = c(1, 0.05,1)
)

pdf(file = paste0(outfolder,"S_sim_sp_i_m_v1.pdf"), width=10, height=8)
S
dev.off()

pdf(file = paste0(outfolder,"S_sim_sp_i_sd_v1.pdf"), width=10, height=3)
row3
dev.off()

# RMSE of CIF and spatial slopes
myscale <- list(bl=0,bu=0.09,bby=0.03,ll=0,uu=0.09)
pA <- to_get_plot_v5(get_plot1_v6,"A",rt=1,GP_label,myscale=myscale,
                     list(GP_cif_m_rmse[,cif_i,],HSGP_cif_m_rmse[,cif_i,]), legend="CIF at Year 5 RMSE")
pB <- to_get_plot_v5(get_plot1_v6,"B",rt=1,GP_label,
                     list(GP_theta1_m_rmse,HSGP_theta1_m_rmse), legend="Spatial slope RMSE")

pC <- to_get_plot_v5(get_plot1_v6,"C",rt=2,GP_label,
                     list(GP_cif_m_rmse[,cif_i,],HSGP_cif_m_rmse[,cif_i,]), legend="CIF at Year 5 RMSE")
pD <- to_get_plot_v5(get_plot1_v6,"D",rt=2,GP_label,
                     list(GP_theta1_m_rmse,HSGP_theta1_m_rmse), legend="Spatial slope RMSE")

S <- plot_grid(
  pA,pB,pC,pD,
  nrow=2, ncol=2,
  rel_widths = c(1,1),
  rel_heights = c(1,1)
)

pdf(file = paste0(outfolder,"S_sim_sp_rmse_v2.pdf"), width=10, height=6)
S
dev.off()

# Posterior coverage

pA <- to_get_plot_v5(get_plot2_v3,"A",rt=1,GP_label,
                     list(GP_cif_coverage[,cif_i,],HSGP_cif_coverage[,cif_i,]), legend="CIF at Year 5 posterior coverage")
pB <- to_get_plot_v5(get_plot2_v3,"B",rt=1,GP_label,
                     list(GP_theta1_coverage,HSGP_theta1_coverage), legend="Spatial slope posterior coverage")
pC <- to_get_plot_v5(get_plot2_v3,"C",rt=2,GP_label,
                     list(GP_cif_coverage[,cif_i,],HSGP_cif_coverage[,cif_i,]), legend="CIF at Year 5 posterior coverage")
pD <- to_get_plot_v5(get_plot2_v3,"D",rt=2,GP_label,
                     list(GP_theta1_coverage,HSGP_theta1_coverage), legend="Spatial slope posterior coverage")

S <- plot_grid(
  pA,pB,pC,pD,
  nrow=2, ncol=2,
  rel_widths = c(1,1),
  rel_heights = c(1,1)
)

pdf(file = paste0(outfolder,"S_sim_sp_cov_v2.pdf"), width=10, height=6)
S
dev.off()

# Baseline hazard rates

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

S <- plot_grid(
  p_GP, p_HSGP,
  labels = c("A", "B"),
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1)
)

pdf(file = paste0(outfolder,"S_sim_bhr_r",j,"_v1.pdf"), width=8, height=3)
S
dev.off()

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

S <- plot_grid(
  p_GP, p_HSGP,
  labels = c("A", "B"),
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1)
)

pdf(file = paste0(outfolder,"S_sim_bhr_r",j,"_v1.pdf"), width=8, height=3)
S
dev.off()

# per simset RMSE of baseline hazard rates

# ylim <- c(0,2)
# lambda_label <- matrix(nrow=k,ncol=m)
# for (j in 1:m) {
#   lambda_label[, j] <- sprintf('lambda["%d,%d"]', j, 1:k)
# }
# lambdabreaks <- c(1,seq(5,k,by=5))
# 
# j <- 1
# pA <- to_get_plot3_v2(list(t(GP_lambda_rmse[,j,]),t(HSGP_lambda_rmse[,j,])),
#                       lambda_label[,j],GP_label,paste0("Risk type ",j),
#                       ylim,xbreaks=lambdabreaks,ylab="RMSE")
# j <- 2 
# pB <- to_get_plot3_v2(list(t(GP_lambda_rmse[,j,]),t(HSGP_lambda_rmse[,j,])),
#                       lambda_label[,j],GP_label,paste0("Risk type ",j),
#                       ylim,xbreaks=lambdabreaks,ylab="RMSE")
# 
# S <- plot_grid(
#   pA, pB,
#   labels = c("A", "B"),
#   label_fontface = "bold",
#   nrow = 2
# )
# 
# pdf(file = paste0(outfolder,"S_sim_bhr_rmse_v1.pdf"), width=12, height=9)
# S
# dev.off()

# Regression coefficients

data.frame(rbind(coxph=coxph1_m_rmse[1:(p+1)],
                 "coxph+groups"=coxph2_m_rmse[1:(p+1)],
                 GP=GP_para_m_rmse[c(1:p,(2*p)+1)],
                 HSGP=HSGP_para_m_rmse[c(1:p,(2*p)+1)])) %>%
  `colnames<-`(paras[c(1:p,(2*p)+1)]) %>%
  kable(digits=2, format="latex")

data.frame(rbind(coxph=coxph1_m_rmse[(p+2):(2*p+2)],
                 "coxph+groups"=coxph2_m_rmse[(p+2):(2*p+2)],
                 GP=GP_para_m_rmse[c((p+1):(2*p),2*p+2)],
                 HSGP=HSGP_para_m_rmse[c((p+1):(2*p),2*p+2)])) %>%
  `colnames<-`(paras[c((p+1):(2*p),2*p+2)]) %>%
  kable(digits=2, format="latex")

# ESS and run time

p1 <- get_plot3(list(GP_ESS,HSGP_ESS),paras,GP_label,"ESS")

RT_df <- data.frame(cbind(value=c(GP_RT,HSGP_RT),
                          label=rep(GP_label, each=nsim) ))

round(c(mean(GP_RT), quantile(GP_RT,c(0.025,0.975)))/60/60,2)
round(c(mean(HSGP_RT), quantile(HSGP_RT,c(0.025,0.975)))/60/60,2)

p2 <- RT_df %>%
  mutate(value=as.numeric(value)) %>%
  ggplot(aes(x=factor(label), y=value, fill=label)) + 
  geom_boxplot() + theme_bw() +
  labs(y="seconds", title="run time", x="") + 
  theme(legend.position="none")

ggarrange(p1,p2,ncol=2,widths=c(6,1))
