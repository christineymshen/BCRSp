# base, no spatial
specnum <- 11

library(tidyverse)

folder <- "application/"

source(paste0(folder,"R/functions.R"))

outcome <- "HP"

df <- readRDS("data/synthetic_v2.rds") %>%
  rename(Longitude=x,Latitude=y)

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

spec <- list(a0=1,b0=1,a1=2,b1=4*10,m=m,n=n,p=dim(X)[2],k=k,X=X,
             s=s,delta=delta,isModeli=isModeli,isModels=isModels,
             data=df_run,outcome=outcome)

if (isModels) spec$W <- W

if (!dir.exists(paste0(folder,"spec"))) dir.create(paste0(folder,"spec"), recursive = TRUE)
saveRDS(spec, paste0(folder,"spec/spec",specnum,".rds"))
