
# update the output folder path here
outfolder <- ""

source(paste0(outfolder, "functions.R"))
nu <- 3/2
m <- 2 # number of risk types
rt_labels <- sprintf("risktype%d",c(1:m))
GP_labels <- c("intercept","slope")

specnum <- 1

isModelrt <- c(T,T) # whether model spatial processes for each risk type
isModeli <- T # whether model spatial intercept, assuming one decision for both risk types
isModels <- T # whether model spatial slope, assuming one decision for both risk types
isJointwithin <- F # whether the spatial intercept and slope within each risk type are jointly modelled
isLMC <- T # whether jointly modeled across risk types, the two entries are respectively for intercept, and slope

# number of chains, number of iterations per chain, number of warm-up iterations, and thinning parameter
nchain <- 4; niter <- 5000; nburn <- 1000; nthin <- 5
pc <- 0.4 # percentage censored in the simulated dataset
kappa0 <- 1; a <- 2; b <- 40; # hyperparameters for the baseline hazard rates
pl<- F # whether the true baseline hazard rates are set to be piecewise linear

bdd <- 1; by <- 0.1 # parameters for the grid for krigging
nGP <- sum(isModeli+isModels) # number of GP per risk type

# get covariates
data_synthetic <- readRDS("synthetic.rds")

d <- cbind(data_synthetic$x, data_synthetic$y)

x1s <- x2s <- seq(-bdd,bdd,by=by)
d_pred <- as.matrix(expand.grid(x1s,x2s)) # kriging locations

# get design matrix and remove intercept
X <- model.matrix(~race+sex+smoking+iswpartner+insurance+age+CMR_readm,data=data_synthetic)[,-1]

# scale and center continuous covariates
idx_cont <- which(apply(X,2,function(x) {!all(x %in% 0:1)}))
X[,idx_cont] <- apply(X[,idx_cont],2,scale)

# center binary covariates
X[,-idx_cont] <- apply(X[,-idx_cont],2,scale,scale=F)

# separate W and X if model spatial slopes for W
if (isModels) {
  idx <- which(colnames(X)=="CMR_readm")
  W <- X[,idx]; X <- X[,-idx]
}

p <- dim(X)[2]
n <- dim(d)[1]; n_pred <- dim(d_pred)[1]
d_all <- rbind(d,d_pred)

set.seed(3)
beta <- matrix(rnorm(p*m,sd=0.5),nrow=p)
beta_w <- c(0.5,-0.7)

# LMC
nw <- 8
nus <- c(1/2,1/2,1/2,3/2,3/2,3/2,5/2,5/2)
ls <- runif(nw,1,2)

d_all_dist <- dist(d_all)
d_all_cov_vec <- vapply(c(1:nw), function(w) vec_matern_cov(d_all_dist,1,nus[w],ls[w]), numeric(length(d_all_dist)))
Ks <- vapply(c(1:nw), function(w) distvec_to_mat(d_all_cov_vec[,w],n+n_pred), matrix(0,n+n_pred,n+n_pred))
LKs <- vapply(c(1:nw), function(w) chol(Ks[,,w]), matrix(0,n+n_pred,n+n_pred))

set.seed(4)
zs <- matrix(rnorm((n+n_pred)*nw),nrow=n+n_pred,ncol=nw)
A <- matrix(runif(4*8,-0.5,0.5),nrow=4)

# A[1,] risk type 1, intercept
# A[2,] risk type 2, intercept
# A[3,] risk type 1, slope
# A[4,] risk type 2, slope

alpha <- matrix(diag(tcrossprod(A)),m,nGP)

ws <- vapply(c(1:nw), function(w) zs[,w] %*% LKs[,,w], numeric(n+n_pred))

# idx for spatial slope
if (isModels) s_idx <- ifelse(isModeli,2,1)
GP_labels_spec <- GP_labels[c(isModeli,isModels)]

f_true <- array(tcrossprod(A,ws),dim=c(m,nGP,n+n_pred),dimnames=list(rt_labels,GP_labels_spec))

# this is used for checking
f_true_check <- array(0,dim=c(m,nGP,n_pred),dimnames=list(rt_labels,GP_labels_spec))
f_true_check <- f_true[,,(n+1):(n+n_pred),drop=F]
if (isModeli){
  f0_level <- rowMeans(f_true_check[,1,])
  f_true_check[,1,] <- f_true_check[,1,] - f0_level
}
if (isModels) f_true_check[,s_idx,] <- f_true_check[,s_idx,] + beta_w

# for baseline hazard rates
gammas <- c(1,1) # scale
alphas <- c(1/0.3,1/0.5) # shape
shift <- c(-10,-5)
k <- 30

# these functions will be used later in checking
lambda_f_sim <- function(t,j,spec){
  spec$gammas[j]*spec$alphas[j]*t^(spec$alphas[j]-1)
}

lambda_f_check <- function(t,j,spec){
  
  t_xbeta_level <- colMeans(spec$X %*% spec$beta)
  spec$lambda_f_sim(t,j,spec)*exp(spec$f0_level[j]+t_xbeta_level[j]+spec$shift[j])
  
}

invert_cdf <- function(U,j,spec){
  
  (U/spec$gammas[j])^(1/spec$alphas[j])
  
}

runspec <- list(isModelrt=isModelrt,isModeli=isModeli,isModels=isModels,isJointwithin=isJointwithin,isLMC=isLMC,
                X=X,W=W,d=d,d_pred=d_pred,n=n,n_pred=n_pred,p=p,k=k,m=m,bdd=1,by=0.1,pc=pc,
                kappa0=kappa0,a=a,b=b,nchain=nchain,niter=niter,nburn=nburn,nthin=nthin,
                nu=nu,alpha=alpha,f_true=f_true,f_true_check=f_true_check,beta=beta,
                pl=pl,lambda_f_sim=lambda_f_sim,lambda_f_check=lambda_f_check,invert_cdf=invert_cdf,
                shift=shift,gammas=gammas,alphas=alphas)

if (isJointwithin) runspec$rho <- rho
if (isModeli) runspec$f0_level <- f0_level
if (isModels) {
  runspec$beta_w <- beta_w
  runspec$s_idx <- s_idx
}

saveRDS(runspec, paste0(outfolder,"spec_",specnum,".rds"))
