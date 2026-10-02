specnum <- 2
runnum <- 1

library(tidyverse)
library(sf)
library(RANN) # clustering
library(cowplot)
library(bayesplot) # betea trace plot
library(kableExtra) # beta tables
library(loo) # cross validation
library(scales) # label_number function

folder <- "application/"

source(paste0(folder,"R/functions.R"))
outfolder <- paste0(folder,"fig/")

datafolder <- "data/"
durham_cities_geo <- readRDS(paste0(datafolder, "GEO/durham_cities_geo.rds"))

load(paste0(datafolder,"GEO/NC_County.RData"))

# rounding scalar
rs <- 10
# number of intervals
ni <- 4

## read in results files
spec <- readRDS(paste0(folder, "spec/spec",specnum,".rds"))

map_sf <- NC_C_2023 %>%
  filter(substr(GEOID,1,5) %in% spec$fips)

if (length(spec$fips)==1 && spec$fips=="37063"){
  durham_cities <- durham_cities_geo %>%
    filter(city=="Durham")
} else {
  durham_cities <- durham_cities_geo %>%
    filter(city != "Carrboro")   
}

res <- readRDS(paste0(folder,"res/spec",specnum,"_run",runnum,".rds"))
fit <- res$res
n_psim <- (res$runspec$niter-res$runspec$nburn)*res$runspec$nchain/res$runspec$nthin

m <- spec$m; n_pred <- spec$n_pred
nGP <- sum(spec$isModeli+spec$isModels)

# extract posterior samples
theta_pred <- array(rstan::extract(fit,"theta_pred")$theta_pred,
                    dim=c(n_psim,m,n_pred,nGP))
lambda <- aperm(rstan::extract(fit, pars="lambda")$lambda, c(1,3,2))

theta0_pred <- theta_pred[,,,1]
p_beta <- rstan::extract(fit,pars="beta")$"beta"

if (spec$isModels){
  p_beta_w <- rstan::extract(fit,pars="beta_w")$"beta_w"
  theta1_pred <- vapply(c(1:n_pred), function(j) theta_pred[,,j,nGP]+p_beta_w, matrix(0,n_psim,m)) #npsim x m x n_pred
  theta1_pred_hr <- exp(theta1_pred)
  theta1_sd <- theta1_m <- matrix(0,n_pred,m)
  
  for (j in 1:m){
    theta1_m[,j] <- colMeans(theta1_pred_hr[,j,])
    theta1_sd[,j] <- apply(theta1_pred_hr[,j,],2,sd)
  }
  
}

### CIF and spatial slopes----

year <- 1
xs <- 52*year

iscompute <- T
if (iscompute){

  # for risk type 1 only
  library(doMC)
  registerDoMC(8)
  
  xbeta_cif <- apply(p_beta,3, function(x) x %*% spec$x_cif)
  theta0_pred_cif <- sweep(theta0_pred,c(1,2),xbeta_cif,"+")
  
  # codes to get summary for CIF
  res_CIF <- foreach(i = 1:n_pred) %dopar% {
    p_lambda_sp_i <- sweep(lambda, MARGIN = c(1, 3), STATS = exp(theta0_pred_cif[,,i]), FUN = "*")
    # n_psim x nxs x m x (p+1)
    res_tmp <- getCIF_v3(xs,spec$s,spec$k,m,p_lambda_sp_i)

    return(res_tmp[,1,,1])
  }
  
  # col is ordered by beta, then by xs
  saveRDS(res_CIF,paste0(folder,"res/CIF_spec", specnum,"run",runnum,"_yr",year,"_v2.rds"))
}

res_CIF <- readRDS(paste0(folder,"res/CIF_spec", specnum,"run",runnum,"_yr",year,"_v2.rds"))

res_CIF_r1 <- sapply(res_CIF, `[`, , 1)
res_CIF_r2 <- sapply(res_CIF, `[`, , 2)

res_CIF_m <- res_CIF_sd <- matrix(0,n_pred,m)

res_CIF_m[,1] <- colMeans(res_CIF_r1)
res_CIF_m[,2] <- colMeans(res_CIF_r2)
res_CIF_sd[,1] <- apply(res_CIF_r1,2,sd)
res_CIF_sd[,2] <- apply(res_CIF_r2,2,sd)

# for risk type 1
j <- 1
m_scale <- list(ll=0,ul=0.14,lby=0.04)
pA <- get_plot1_v11(spec$geo_pred_sf,map_sf,res_CIF_m[,j],center=mean(res_CIF_m[,j]),
                    legend="CIF at Year 1 posterior mean",m_scale=m_scale,text_df=durham_cities)

m_scale <- get_plot_scale_v5(rs,ni,list(res_CIF_sd[,j]))
pB <- get_plot1_v11(spec$geo_pred_sf,map_sf,res_CIF_sd[,j],
                    legend="CIF at Year 1 posterior SD",m_scale=m_scale,text_df=durham_cities)

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

# figure for the draft paper
pdf(file = paste0(outfolder,"paper/spec",specnum,"_run",runnum,"_S_sp_r",j,"_v2.pdf"), width=10, height=10)
Final
dev.off()


### Clustering ----

# cluster CIFs

k <- 3
labels <- c("Low","Medium","High")

j <- 1
cluster_k <- kmeans(res_CIF_m[,j],k)
df <- data.frame(cbind(cluster=cluster_k$cluster,
                       m = res_CIF_m[,j]))

df_label <- df %>%
  group_by(cluster) %>%
  summarise(m_mean=round(mean(m),2)) %>%
  arrange(m_mean) %>%
  mutate(label=labels)

df_label

df_i <- df %>%
  left_join(df_label, join_by("cluster"))

# computes posterior uncertainty given the centers
i_cluster_pred <- vapply(c(1:n_psim), function(i) nn2(cluster_k$centers,res_CIF_r1[i,],k=1)$nn.idx,
                         integer(n_pred))
i_cluster_prop <- rowMeans(i_cluster_pred == cluster_k$cluster)


pA <- print(get_plot3_v1(spec$geo_pred_sf, map_sf, factor(df_i$label, levels=labels),
             legend="CIF at Year 1 clusters",
             text_df = durham_cities,
             legendsize=8,textsize=3))

m_scale <- list(ll=0,ul=1,lby=0.2)
pB <- print(get_plot1_v11(spec$geo_pred_sf, map_sf, i_cluster_prop, m_scale=m_scale,
              legend="CIF at Year 1 cluster posterior probability",
              text_df = durham_cities))

# cluster spatial slopes
cluster_k <- kmeans(theta1_m[,j],k)

df <- data.frame(cbind(cluster=cluster_k$cluster,
                       m = theta1_m[,j]))

df_label <- df %>%
  group_by(cluster) %>%
  summarise(m_mean=round(mean(m),2)) %>%
  arrange(m_mean) %>%
  mutate(label=labels)

df_label

df_i <- df %>%
  left_join(df_label, join_by("cluster"))

# computes posterior uncertainty given the centers
i_cluster_pred <- vapply(c(1:n_psim), function(i) nn2(cluster_k$centers,theta1_pred_hr[i,1,],k=1)$nn.idx,integer(n_pred))
i_cluster_prop <- rowMeans(i_cluster_pred == cluster_k$cluster)

pC <- print(get_plot3_v1(spec$geo_pred_sf, map_sf, factor(df_i$label, levels=labels),
                         legend="Spatial slope clusters",
                         text_df = durham_cities,
                         legendsize=8,textsize=3))

m_scale <- list(ll=0,ul=1,lby=0.2)
pD <- print(get_plot1_v11(spec$geo_pred_sf, map_sf, i_cluster_prop, m_scale=m_scale,
                          legend="Slope cluster posterior probability",
                          text_df = durham_cities))

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
pdf(file = paste0(outfolder,"paper/spec",specnum,"_run",runnum,"_S_clust_r",j,"_v2.pdf"), width=10, height=10)
Final
dev.off()


### Beta ----

# extract posterior samples
p_beta <- as.matrix(fit,pars="beta")
rowspeclabels <- c("Race - Black", "Race - Other", "Sex - Female", "Smoking - Former", 
                   "Smoking - Current", "No partner", "Insurance - Medicaid",
                   "Insurance - Medicare", "Insurance - Other", "Age")
risklabels <- c("Readmission risk", "Mortality risk")

## regression coefficients plots

par(mfrow=c(1,2),mar=c(3,1,2,1), mgp=c(1.75, 0.75, 0), oma=c(1,10,1,1))
HR <- vector("list", length=spec$m)
for (j in 1:m){
  if (j==1){
    HR[[j]] <- getHR_multi_v3(p_beta,rowspeclabels,risktype=j)
  } else {
    HR[[j]] <- getHR_multi_v3(p_beta,rowspeclabels,risktype=j, ordered_labels = HR[[1]]$label)
  }
  values <- HR[[j]][,3]
  dfs <- list(HR[[j]])
  xmax <- ceiling(max(values)*rs)/rs
  
  if (j==1){
    plotHRs_v3(dfs, risklabels[j], cex=1.2, xlim=c(0,xmax))
  } else {
    plotHRs_v3(dfs, risklabels[j], yaxt=F, cex=1.2, xlim=c(0,xmax))
  }
}


### Baseline hazard rates ----

# n_psim x n x m
theta0 <- rstan::extract(fit,"theta0")$theta0
theta0_levels <- apply(theta0, MARGIN = c(1,3), mean)

j <- 1
lambda[,,j] <- lambda[,,j] * exp(theta0_levels[,j])

lambda_m <- colMeans(lambda[,,j])
lambda_u <- apply(lambda[,,j],2,quantile,0.975)
lambda_l <- apply(lambda[,,j],2,quantile,0.025)
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
lambda[,,j] <- lambda[,,j] * exp(theta0_levels[,j])
lambda_m <- colMeans(lambda[,,j])
lambda_u <- apply(lambda[,,j],2,quantile,0.975)
lambda_l <- apply(lambda[,,j],2,quantile,0.025)
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
pdf(file = paste0(outfolder,"paper/spec",specnum,"_run",runnum,"_S_bhr_v1.pdf"), width=10.5, height=3.5)
Final
dev.off()


### Spatial intercepts ----

for (j in 1:m){
  theta0_pred[,j,] <- theta0_pred[,j,] - theta0_levels[,j]
}

theta0_m <- theta0_sd <- matrix(0,n_pred,m)

for (j in 1:m){
  theta0_m[,j] <- colMeans(theta0_pred[,j,])
  theta0_sd[,j] <- apply(theta0_pred[,j,],2,sd)
}

# note that pdf codes don't work if I put them into a for loop
## risk type 1
j <- 1
m_scale <- list(ll=-0.6,ul=0.7,lby=0.3)
pA <- print(get_plot1_v11(spec$geo_pred_sf,map_sf,theta0_m[,j], center=0,
                          legend="Readmission risk intercept posterior mean",m_scale=m_scale,text_df=durham_cities))

m_scale <- get_plot_scale_v2(rs,ni,j,list(theta0_sd))
pB <- print(get_plot1_v11(spec$geo_pred_sf,map_sf,theta0_sd[,j],
                          legend="Readmission risk intercept posterior SD",m_scale=m_scale,text_df=durham_cities))

j <- 2
m_scale <- get_plot_scale_v2(rs,ni,j,list(theta0_m))
m_scale <- list(ll=-0.1,ul=0.15,lby=0.05)
pC <- print(get_plot1_v11(spec$geo_pred_sf,map_sf,theta0_m[,j],center=0,
                          legend="Mortality risk intercept posterior mean",m_scale=m_scale,text_df=durham_cities))


m_scale <- get_plot_scale_v2(rs,ni,j,list(theta0_sd))
pD <- print(get_plot1_v11(spec$geo_pred_sf,map_sf,theta0_sd[,j],
                          legend="Mortality risk intercept posterior SD",m_scale=m_scale,text_df=durham_cities))


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

# figure for the draft paper
pdf(file = paste0(outfolder,"paper/spec",specnum,"_run",runnum,"_S_sp_i_v1.pdf"), width=10, height=10)
Final
dev.off()

