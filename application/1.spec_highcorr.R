library(tidyverse)
library(sp)
library(sf)

wfolder <- ""
folder <- ""
source(paste0(wfolder, "R/functions.R"))

outcomes <- "HP"; specnum <- 3

load(paste0(folder,"Processed/data_20241024.RData"))
load(paste0(folder,"GEO/NC_BG.RData"))

# Chatham, Orange, Person, Granville, Wake, Durham
durham_NH_fips <- c("37037","37135","37145","37077","37183","37063")

df <- data %>%
  filter(substr(NC_FIPS_patch,1,5) %in% durham_NH_fips) %>%
  select(age,race,sex,iswpartner,insurance,CMR_readm, bp, Longitude, Latitude,
         mhi_adj, bdh, hsl, ue, bb, SVI, ADI1, LILATracts_Vehicle, smoking,
         time_c, time_ed, time_d,time_hp,time_fr,time_fa, FIPS_20) %>%
  filter(race != "NA", iswpartner != "NA", smoking != "NA", 
         !ADI1 %in% c("GQ-PH","QDI","PH","GQ","NONE","NA")) %>%
  mutate(age = as.numeric(age),
         ADIc = as.numeric(ADI1),
         ADI_cat = case_when(ADIc>= 1 & ADIc<=15 ~ "low",
                             ADIc>=16 & ADIc<=85 ~ "medium",
                             ADIc>=86 & ADIc<=100 ~ "high"), 
         SVI_cat = case_when(SVI<=0.25 ~ "low",
                             SVI>0.25 & SVI<=0.5 ~ "low-medium",
                             SVI>0.5 & SVI<=0.75 ~ "medium",
                             SVI>0.75 & SVI<=1 ~ "high"),
         race = factor(race,levels=c("white","black","other")),
         sex = factor(sex,levels=c("Male","Female")),
         iswpartner = factor(iswpartner, levels=c("yes","no")),
         insurance = factor(insurance, levels=c("Government","Commercial","WCSC","Self-Pay")),
         smoking = factor(smoking, levels=c("Never","Former","Smoker")),
         ADI_cat=factor(ADI_cat, levels=c("low","medium","high")),
         SVI_cat=factor(SVI_cat, levels=c("low", "low-medium", "medium","high"))) %>%
  rename(SVIc=SVI,LILA=LILATracts_Vehicle,fips=FIPS_20)

# UTM projection

# get geo coordinates
geo_raw <- geo <- df %>%
  select(Longitude,Latitude) %>%
  filter(if_any(everything(), ~ !is.na(.)))

# try to converge lat and long to x and y
coordinates(geo) <- c("Longitude","Latitude")

# actually not sure what EPSG codes for CRS, WGS84 is said to be the most commonly used one
proj4string(geo) <- CRS("+proj=longlat +datum=WGS84")
geo_proj <- as.data.frame(spTransform(geo, CRS("+proj=utm +zone=17 ellps=WGS84"))) %>%
  `colnames<-`(c("x","y"))

# euclidean distances
dist_geo_proj <- as.matrix(dist(geo_proj, diag=T,upper=T))

# vector version
dist_geo_proj_v <- c(dist_geo_proj[lower.tri(dist_geo_proj)])

# compute distances based on lat and long and check
# maximum difference is 335m, ok
# mat <- distm(geo_raw, geo_raw, fun = distHaversine)
# mat_v <- c(mat[lower.tri(mat)])
# max(abs(mat_v-dist_geo_proj_v))

# center the x and y's for HSGP use
x_adj <- round((max(geo_proj$x)-min(geo_proj$x))/2,0) + min(geo_proj$x)
y_adj <- round((max(geo_proj$y)-min(geo_proj$y))/2,0) + min(geo_proj$y)
scale_x <- 10000; scale_y <- 10000

d <- geo_proj %>%
  mutate(x=(x-x_adj)/scale_x,
         y=(y-y_adj)/scale_y)

dist_d <- as.matrix(dist(d, diag=T,upper=T))
dist_d_v <- c(dist_d[lower.tri(dist_d)])

# replace lat and long with x and y
df <- cbind(df, d) %>%
  select(-c(Longitude, Latitude))

# range(df_run$x); range(df_run$y)

S1 <- 5.4; S2 <- 5.6

xs <- seq(-S1,S1,by=0.3)
ys <- seq(-S2,S2,by=0.3)

d_pred <- expand.grid(xs,ys)

geo_pred_proj <- d_pred %>%
  `colnames<-`(c("x","y")) %>%
  mutate(x = x*scale_x + x_adj,
         y = y*scale_y + y_adj)

coordinates(geo_pred_proj) <- ~x+y
proj4string(geo_pred_proj) <- CRS("+proj=utm +zone=17 ellps=WGS84")

geo_pred <- spTransform(geo_pred_proj, CRS("+proj=longlat +datum=WGS84"))
geo_pred <- as.data.frame(geo_pred) %>%
  `colnames<-`(c("long", "lat"))

durham_NH_BG <- NC_BG_2020 %>%
  filter(substr(GEOID,1,5) %in% durham_NH_fips)

geo_pred_sf <- st_as_sf(geo_pred, coords = c("long", "lat"), crs = st_crs(durham_NH_BG))
geo_sf <- st_as_sf(geo, coords = c("long", "lat"), crs = st_crs(durham_NH_BG))

in_BG <- st_within(geo_pred_sf, durham_NH_BG, sparse=FALSE)

is_in_durham <- apply(in_BG,1,any)
d_pred <- d_pred[is_in_durham,]
geo_pred_sf <- geo_pred_sf[is_in_durham,]
is_in_BG <- in_BG[is_in_durham,]

idx <- apply(is_in_BG,1,which)
pred_fips <- durham_NH_BG$GEOID[idx]

isModels <- T

df_run <- get_df_run(df,outcome)

X <- model.matrix(~race+sex+smoking+iswpartner+insurance+age+CMR_readm,data=df_run)[,-1]
idx_scale <- which(apply(X,2,function(x) {!all(x %in% 0:1)}))
X[,idx_scale] <- apply(X[,idx_scale],2,scale)
X[,-idx_scale] <- apply(X[,-idx_scale],2,scale,scale=F)

if (isModels){
  idx <- which(colnames(X)=="CMR_readm")
  W <- X[,idx]; X <- X[,-idx]
}

max_time <- ceiling(max(df_run$time))
k <- 50
s <- seq(0,max_time,length.out=k+1)
n <- dim(df_run)[1]; m <- 2; n_pred <- dim(d_pred)[1]
delta <- matrix(0,n,m)
idx <- cbind(c(1:n), df_run$event)[df_run$event!=0,]
delta[idx] <- 1

spec <- list(a0=0.1,b0=0.1,a1=2,b1=80,tau=4,al=2,bl=1,m=m,n=n,n_pred=n_pred,p=dim(X)[2],k=k,X=X,
             s=s,delta=delta,d=d,d_pred=d_pred,isModeli=T,isModels=isModels,
             data=df_run,S=c(S1,S2),geo_pred_sf=geo_pred_sf,geo_sf=geo_sf,
             outcome=outcome,fips=durham_NH_fips,pred_fips=pred_fips)

if (isModels) spec$W <- W

saveRDS(spec, paste0(wfolder,"spec/spec",specnum,".rds"))
