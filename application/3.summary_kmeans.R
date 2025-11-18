
specnums <- 163
runnums <- 1

library(tidyverse)
library(sf)
library(viridis)
library(RANN)
library(cowplot)

outfolder <- ""

wdir <- ""

folder <- paste0(wdir, "CRS/")
source(paste0(folder,"R/functions.R"))

datafolder <- paste0(wdir, "Data/")
durham_cities_geo <- readRDS(paste0(datafolder, "GEO/durham_cities_geo.rds"))

load(paste0(datafolder,"GEO/NC_CT.RData"))
load(paste0(datafolder,"GEO/NC_County.RData"))

risklabels <- c("readmission risk", "mortality risk")

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

for (j in 1:m){
  theta0_pred_levels[,j] <- rowMeans(theta0_pred[,j,])
  theta0_pred[,j,] <- theta0_pred[,j,] - theta0_pred_levels[,j]
  theta0_m[,j] <- colMeans(theta0_pred[,j,])
  theta0_sd[,j] <- apply(theta0_pred[,j,],2,sd)
}

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

nn.idx <- nn2(as.matrix(spec$d_pred), as.matrix(spec$d),k=1)$nn.idx

# spatial intercepts

k <- 3
cluster_k <- kmeans(theta0_m[,1],k)

labels <- c("Low","Medium","High")

df <- data.frame(cbind(cluster=cluster_k$cluster,
                       m = theta0_m[,1]))

df_label <- df %>%
  group_by(cluster) %>%
  summarise(m_mean=round(mean(m),2)) %>%
  arrange(m_mean) %>%
  mutate(label=labels)

df_i <- df %>%
  left_join(df_label, join_by("cluster"))

# computes posterior uncertainty given the centers
i_cluster_pred <- vapply(c(1:n_psim), function(i) nn2(cluster_k$centers,theta0_pred[i,1,],k=1)$nn.idx,
                       integer(n_pred))
i_cluster_prop <- rowMeans(i_cluster_pred == cluster_k$cluster)

pA <- print(ggplot() +
  geom_sf(data=spec$geo_pred_sf, aes(color=factor(df_i$label, levels=labels)),
          shape=15, size=10) +
  geom_sf(data=map_sf, fill=NA, color="#bdbdbd", size=1) +
  theme_bw() +
  scale_color_viridis_d(option="A", direction=-1, begin=0.6) +
  theme(panel.border = element_blank(),
        legend.position = "bottom",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.margin = margin(t = -15, unit = "pt")) +
  geom_text(data = durham_cities, aes(lon, lat, label = city), size = 3, fontface="bold")+
  labs(x="",y="",color="Spatial intercept clusters"))

pB <- ggplot() +
  geom_sf(data=spec$geo_pred_sf, aes(color=i_cluster_prop),
          shape=15, size=10) +
  geom_sf(data=map_sf, fill=NA, color="#bdbdbd", size=1) +
  theme_bw() +
  scale_color_viridis_c(option="A", direction=-1, begin=0.5,
                        limits=c(0,1), breaks=seq(0,1,by=0.2)) +
  theme(panel.border = element_blank(),
        legend.position = "bottom",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.margin = margin(t = -15, unit = "pt")) +
  geom_text(data = durham_cities, aes(lon, lat, label = city), size = 3, fontface="bold")+
  labs(x="",y="",color="Intercept cluster posterior probability")

# spatial slopes

k <- 3
cluster_k <- kmeans(theta1_m[,1],k)

df <- data.frame(cbind(cluster=cluster_k$cluster,
                       m = theta1_m[,1]))

df_label <- df %>%
  group_by(cluster) %>%
  summarise(m_mean=round(mean(m),2)) %>%
  arrange(m_mean) %>%
  mutate(label=labels)

df_s <- df %>%
  left_join(df_label, join_by("cluster"))

# computes posterior uncertainty given the centers
s_cluster_pred <- vapply(c(1:n_psim), function(i) nn2(cluster_k$centers,theta1_pred_hr[i,1,],k=1)$nn.idx,
                       integer(n_pred))
s_cluster_prop <- rowMeans(s_cluster_pred == cluster_k$cluster)

pC <- print(ggplot()+
  geom_sf(data=spec$geo_pred_sf, aes(color=factor(df_s$label, levels=labels)),
          shape=15, size=10) +
  geom_sf(data=map_sf, fill=NA, color="#bdbdbd", size=1) +
  theme_bw() +
  scale_color_viridis_d(option="A", direction=-1, begin=0.6) +
  theme(panel.border = element_blank(),
        legend.position = "bottom",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.margin = margin(t = -15, unit = "pt")) +
  geom_text(data = durham_cities, aes(lon, lat, label = city), size = 3, fontface="bold") +
  labs(x="",y="",color="Spatial slope clusters"))

pD <- print(ggplot()+
  geom_sf(data=spec$geo_pred_sf, aes(color=s_cluster_prop),
          shape=15, size=10) +
  geom_sf(data=map_sf, fill=NA, color="#bdbdbd", size=1) +
  theme_bw() +
  scale_color_viridis_c(option="A", direction=-1, begin=0.4,
                        limits=c(0,1), breaks = seq(0, 1, by = 0.2)) +
  theme(panel.border = element_blank(),
        legend.position = "bottom",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.margin = margin(t = -15, unit = "pt")) +
  geom_text(data = durham_cities, aes(lon, lat, label = city), size = 3, fontface="bold") +
  labs(x="",y="",color="Slope cluster posterior probability"))

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

pdf(file = paste0(outfolder,"/fig/F9_clustering_r1_v4.pdf"), width=10, height=10)
Final
dev.off()

