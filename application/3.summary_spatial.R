specnums <- 163
runnums <- 1

library(tidyverse)
library(sf)
library(viridis)
library(RANN)
library(cowplot)
library(scales)

outfolder <- ""

wdir <- ""

folder <- paste0(wdir, "CRS/")
source(paste0(folder,"R/functions.R"))

datafolder <- paste0(wdir, "Data/")
durham_cities_geo <- readRDS(paste0(datafolder, "GEO/durham_cities_geo.rds"))

load(paste0(datafolder,"GEO/NC_CT.RData"))
load(paste0(datafolder,"GEO/NC_County.RData"))

risklabels <- c("Readmission risk", "Mortality risk")

# rounding scalar
rs <- 10
# number of intervals
ni <- 4
qtiles <- c(0.05,0.95)

## read in results files
spec <- readRDS(paste0(folder, "spec/spec",specnums,".rds"))

# in the past, we used census tract boundaries for plotting
# map_sf <- NC_CT_2020 %>%
#   filter(substr(GEOID,1,5) %in% spec$fips)

map_sf <- NC_C_2023 %>%
  filter(substr(GEOID,1,5) %in% spec$fips)

if (length(spec$fips)==1 && spec$fips=="37063"){
  durham_cities <- durham_cities_geo %>%
    filter(city=="Durham")
} else {
  durham_cities <- durham_cities_geo %>%
    filter(city != "Carrboro")   
}

fit <- readRDS(paste0(folder,"res/HSGP_spec",specnums,"_run",runnums,".rds"))$res

m <- spec$m; n_pred <- spec$n_pred
tmp <- readRDS(paste0(folder,"res/HSGP_spec",specnums,"_run",runnums,".rds"))
n_psim <- (tmp$runspec$niter-tmp$runspec$nburn)*tmp$runspec$nchain/tmp$runspec$nthin

nGP <- sum(spec$isModeli+spec$isModels)
theta_pred <- array(rstan::extract(fit,"theta_pred")$theta_pred,
                    dim=c(n_psim,m,n_pred,nGP))
theta0_pred <- theta_pred[,,,1]
theta0_pred_levels <- matrix(0,n_psim,m)
theta0_m <- theta0_sd <- matrix(0,n_pred,m)
#theta0_u <- theta0_l <- theta0_m_masked <- matrix(0,n_pred,m)

for (j in 1:m){
  theta0_pred_levels[,j] <- rowMeans(theta0_pred[,j,])
  theta0_pred[,j,] <- theta0_pred[,j,] - theta0_pred_levels[,j]
  theta0_m[,j] <- colMeans(theta0_pred[,j,])
  theta0_sd[,j] <- apply(theta0_pred[,j,],2,sd)
  #theta0_l[,j] <- apply(theta0_pred[,j,],2,quantile,qtiles[1])
  #theta0_u[,j] <- apply(theta0_pred[,j,],2,quantile,qtiles[2])
}
#theta0_mask <- theta0_l >0 | theta0_u < 0

if (spec$isModels){
  p_beta_w <- as.matrix(fit,pars="beta_w")
  theta1_pred <- vapply(c(1:n_pred), function(j) theta_pred[,,j,nGP]+p_beta_w, matrix(0,n_psim,m)) #npsim x m x n_pred
  theta1_pred_hr <- exp(theta1_pred)
  theta1_sd <- theta1_m <- theta1_u <- theta1_l <- theta1_m_masked <- matrix(0,n_pred,m)
  
  for (j in 1:m){
    theta1_m[,j] <- colMeans(theta1_pred_hr[,j,])
    theta1_sd[,j] <- apply(theta1_pred_hr[,j,],2,sd)
    theta1_l[,j] <- apply(theta1_pred_hr[,j,],2,quantile,qtiles[1])
    theta1_u[,j] <- apply(theta1_pred_hr[,j,],2,quantile,qtiles[2])
  }
  theta1_mask <- theta1_l >1 | theta1_u < 1
  theta1_m_masked[theta1_mask] <- theta1_m[theta1_mask]
}


## Plots for spatial surfaces

## note that pdf codes don't work if I put them into a for loop

## risk type 1
j <- 1

## Spatial intercepts

m_scale <- get_plot_scale_v2(rs,ni,j,list(theta0_m))
pA <- get_plot1_v11(spec$geo_pred_sf,map_sf,theta0_m[,j],
              legend="Intercept posterior mean",m_scale=m_scale,text_df=durham_cities)


m_scale <- get_plot_scale_v2(rs,ni,j,list(theta0_sd))
pB <- get_plot1_v11(spec$geo_pred_sf,map_sf,theta0_sd[,j],
              legend="Intercept posterior SD",m_scale=m_scale,text_df=durham_cities)


## Spatial slopes
m_scale <- get_plot_scale_v2(rs,ni,j,list(theta1_m))
pC <- get_plot1_v11(spec$geo_pred_sf,map_sf,theta1_m[,j],
                    center=mean(theta1_m[,j]),
              legend="Slope posterior mean",m_scale=m_scale,text_df=durham_cities)


m_scale <- get_plot_scale_v2(rs,ni,j,list(theta1_sd))
pD <- get_plot1_v11(spec$geo_pred_sf,map_sf,theta1_sd[,j], 
              legend="Slope posterior SD",m_scale=m_scale,text_df=durham_cities)


row1 <- plot_grid(
  pA, pB,
  labels = c("A", "B"),
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1)
)
row2 <- plot_grid(
  pC, pD,
  labels = c("C", "D"),
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1)
)

Final <- plot_grid(
  row1,
  row2,
  ncol = 1,
  rel_heights = c(1, 1)
)
pdf(file = paste0(outfolder,"/fig/F8_r",j,"_v3.pdf"), width=10, height=10)
Final
dev.off()

j <- 2

## Spatial intercepts

m_scale <- get_plot_scale_v2(rs,ni,j,list(theta0_m))
pA <- get_plot1_v11(spec$geo_pred_sf,map_sf,theta0_m[,j],
                    legend="Intercept posterior mean",m_scale=m_scale,text_df=durham_cities)


m_scale <- get_plot_scale_v2(rs,ni,j,list(theta0_sd))
pB <- get_plot1_v11(spec$geo_pred_sf,map_sf,theta0_sd[,j],
                    legend="Intercept posterior SD",m_scale=m_scale,text_df=durham_cities)


## Spatial slopes
m_scale <- get_plot_scale_v2(rs,ni,j,list(theta1_m))
pC <- get_plot1_v11(spec$geo_pred_sf,map_sf,theta1_m[,j],
                    center=mean(theta1_m[,j]),
                    legend="Slope posterior mean",m_scale=m_scale,text_df=durham_cities)


m_scale <- get_plot_scale_v2(rs,ni,j,list(theta1_sd))
pD <- get_plot1_v11(spec$geo_pred_sf,map_sf,theta1_sd[,j],
                    legend="Slope posterior SD",m_scale=m_scale,text_df=durham_cities)


row1 <- plot_grid(
  pA, pB,
  labels = c("A", "B"),
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1)
)
row2 <- plot_grid(
  pC, pD,
  labels = c("C", "D"),
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1)
)

Final <- plot_grid(
  row1,
  row2,
  ncol = 1,
  rel_heights = c(1, 1)
)
pdf(file = paste0(outfolder,"/fig/S_r",j,"_v3.pdf"), width=10, height=10)
Final
dev.off()




## Baseline hazard rates

lambda <- rstan::extract(fit, pars="lambda")$lambda
j <- 1
lambda[,j,] <- lambda[,j,] * exp(theta0_pred_levels[,j])

lambda_m <- colMeans(lambda[,j,])
lambda_u <- apply(lambda[,j,],2,quantile,0.975)
lambda_l <- apply(lambda[,j,],2,quantile,0.025)
dfA <- data.frame(m=lambda_m, l=lambda_l, u=lambda_u)

pA <- ggplot(dfA, aes(x = spec$s[-1]/52, y = m)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = l, ymax = u),width = 0.2) +
  labs(x = "Years from initial surgery", y = "Baseline hazard rate", 
       title= "Readmission risk") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())

j <- 2
lambda[,j,] <- lambda[,j,] * exp(theta0_pred_levels[,j])

lambda_m <- colMeans(lambda[,j,])
lambda_u <- apply(lambda[,j,],2,quantile,0.975)
lambda_l <- apply(lambda[,j,],2,quantile,0.025)
dfB <- data.frame(m=lambda_m, l=lambda_l, u=lambda_u)

pB <- ggplot(dfB, aes(x = spec$s[-1]/52, y = m)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = l, ymax = u),width = 0.2) +
  labs(x = "Years from initial surgery", y = "Baseline hazard rate",
       title= "Mortality risk") +
  theme_bw() +
  scale_y_continuous(labels = label_number()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())

Final <- plot_grid(
  pA, pB,
  labels = c("A", "B"),
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1)
)
pdf(file = paste0(outfolder,"/fig/F7_bhr_v4.pdf"), width=10.5, height=3.5)
Final
dev.off()

