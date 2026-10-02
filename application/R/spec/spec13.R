# high correlation, no spatial
specnum <- 13

library(tidyverse)

folder <- "application/"

source(paste0(folder,"R/functions.R"))

outcome <- "HP"

load(paste0(folder,"../Data/Processed/data_20260524.RData"))

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
         insurance = factor(insurance, levels=c("Commercial","Medicaid","Medicare","Other")),
         smoking = factor(smoking, levels=c("Never","Former","Smoker")),
         ADI_cat=factor(ADI_cat, levels=c("low","medium","high")),
         SVI_cat=factor(SVI_cat, levels=c("low", "low-medium", "medium","high"))) %>%
  rename(SVIc=SVI,LILA=LILATracts_Vehicle,fips=FIPS_20)

isModels <- F
isModeli <- F
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
n <- dim(df_run)[1]; m <- 2
delta <- matrix(0,n,m)
idx <- cbind(c(1:n), df_run$event)[df_run$event!=0,]
delta[idx] <- 1

spec <- list(a0=0.1,b0=0.1,a1=2,b1=80,m=m,n=n,p=dim(X)[2],k=k,X=X,
             s=s,delta=delta,isModeli=isModeli,isModels=isModels,
             data=df_run,outcome=outcome,fips=durham_NH_fips)

if (isModels) spec$W <- W

saveRDS(spec, paste0(folder,"spec/spec",specnum,".rds"))
