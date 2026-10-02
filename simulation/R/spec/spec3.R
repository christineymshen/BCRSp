# n=100

specnum <- 3
dataspecnum <- 1

folder <- "simulation/"
source(paste0(folder,"R/functions.R"))
dataspec <- readRDS(paste0("application/spec/spec",dataspecnum,".rds"))

cut_x_l <- -1.52
cut_x_u <- 1.52
cut_y_l <- -1.52
cut_y_u <- 1.52

# used the same scaling as in spec1 so that the map would be on the same scale
data <- dataspec$data %>%
  mutate(num=row_number()) %>%
  filter(x<=cut_x_u,x>=cut_x_l,y<=cut_y_u,y>=cut_y_l) %>%
  mutate(x=x/2.4,y=y/2.4)

# range(data$x)
# range(data$y)

X <- dataspec$X[data$num,]
idx_cont <- which(apply(X,2,function(x) length(unique(x))!=2))
X[,idx_cont] <- apply(X[,idx_cont,drop=F],2,scale)
X[,-idx_cont] <- apply(X[,-idx_cont],2,scale,scale=F)
x_cif <- apply(X,2,min)
x_cif[idx_cont] <- 0
W <- as.numeric(scale(dataspec$W[data$num]))

bdd_x <- 1.6
bdd_y <- 1.6
by_x <- 0.1
by_y <- 0.1
x1s <- seq(-bdd_x,bdd_x,by=by_x)
x2s <- seq(-bdd_y,bdd_y,by=by_y)

d_pred <- expand.grid(x1s,x2s)
d <- data %>%
  select(x,y)

pc <- 0.4 # percentage of censoring
pl <- F
isModeli <- T # whether model spatial intercepts
isModels <- T # whether model spatial slopes
isLMC <- T

nu <- 3/2

nGP <- sum(isModeli+isModels) # number of GP per risk type
m <- dataspec$m # number of risk types
p <- dim(X)[2]
n <- dim(data)[1]; n_pred <- dim(d_pred)[1]
names(d_pred) <- names(d)
d_all <- rbind(d,d_pred)

set.seed(1)
beta <- matrix(rnorm(p*m,sd=0.5),nrow=p)
beta_w <- c(0.5,-0.7)

# LMC
nugget <- 0.0001
nw <- 8
nus <- c(1/2,1/2,1/2,3/2,3/2,3/2,5/2,5/2)
ls <- runif(nw,1,2)

d_all_dist <- dist(d_all)
d_all_cov_vec <- vapply(c(1:nw), function(w) vec_matern_cov(c(d_all_dist),1,nus[w],ls[w]), numeric(length(d_all_dist)))
Ks <- vapply(c(1:nw), function(w) distvec_to_mat(d_all_cov_vec[,w],n+n_pred,nugget=nugget), matrix(0,n+n_pred,n+n_pred))
LKs <- vapply(c(1:nw), function(w) chol(Ks[,,w]), matrix(0,n+n_pred,n+n_pred))

set.seed(2)
zs <- matrix(rnorm((n+n_pred)*nw),nrow=n+n_pred,ncol=nw)
A <- matrix(runif((m+nGP)*nw,-0.5,0.5),nrow=m+nGP)

rt_labels <- sprintf("risktype%d",c(1:m))
GP_labels <- c("intercept","slope")

# A[1,] risk type 1, intercept
# A[2,] risk type 2, intercept
# A[3,] risk type 1, slope
# A[4,] risk type 2, slope

alpha <- sqrt(matrix(diag(tcrossprod(A)),m,nGP))

ws <- vapply(c(1:nw), function(w) zs[,w] %*% LKs[,,w], numeric(n+n_pred))

# idx for spatial slope
if (isModels) s_idx <- ifelse(isModeli,2,1)
GP_labels_spec <- GP_labels[c(isModeli,isModels)]

f_true <- array(tcrossprod(A,ws),dim=c(m,nGP,n+n_pred),dimnames=list(rt_labels,GP_labels_spec))

# this is used for checking
f_true_check <- array(0,dim=c(m,nGP,n_pred),dimnames=list(rt_labels,GP_labels_spec))
f_true_check <- f_true[,,(n+1):(n+n_pred),drop=F]
if (isModeli){
  # f0_level <- rowMeans(f_true_check[,1,])
  f0_level <- rowMeans(f_true[,1,1:n])
  f_true_check[,1,] <- f_true_check[,1,] - f0_level
}
if (isModels) f_true_check[,s_idx,] <- f_true_check[,s_idx,] + beta_w

vis_surface(f_true_check,d_pred,GP_labels_spec)

# for baseline hazard rates
gammas <- c(1,1) # scale
alphas <- c(1/0.3,1/0.5) # shape
shift <- c(-10,-5)
k <- 30

# this function is used to simulate the datasets
# I can't remember why I chose to put spec$shift into xbeta though
lambda_f_sim <- function(t,j,spec){
  spec$gammas[j]*spec$alphas[j]*t^(spec$alphas[j]-1)
}

# this function is used to check the baseline hazard rates through post processing
# if the covariates are not centered, also should include the level of xbeta (see previous versions of spec.R scripts)
# here I removed xbeta because I have centered all the covariates
lambda_f_check <- function(t,j,spec){
  spec$lambda_f_sim(t,j,spec)*exp(spec$f0_level[j]+spec$shift[j])
}

invert_cdf <- function(U,j,spec){
  (U/spec$gammas[j])^(1/spec$alphas[j])
}

spec <- list(isModeli=isModeli,isModels=isModels, isLMC=isLMC,
             X=X,W=W,d=d,d_pred=d_pred,n=n,n_pred=n_pred,p=p,k=k,m=m,bdd_x=bdd_x,bdd_y=bdd_y,
             pc=pc,a0=dataspec$a0,b0=dataspec$b0,a1=dataspec$a1,b1=dataspec$b1,
             al=dataspec$al,bl=dataspec$bl,tau=dataspec$tau,
             nu=nu,alpha=alpha,f_true=f_true,f_true_check=f_true_check,beta=beta,
             pl=pl,lambda_f_sim=lambda_f_sim,lambda_f_check=lambda_f_check,invert_cdf=invert_cdf,
             shift=shift,gammas=gammas,alphas=alphas,ncenters=4,
             nchain=4,niter=5000,nburn=1000,nthin=5)

if (isModeli) spec$f0_level <- f0_level
if (isModels) {
  spec$beta_w <- beta_w
  spec$s_idx <- s_idx
}

# CIF
cif_xs <- c(1:10)
t_cif <- getcif_pc_v2(cif_xs,seq(0,11,by=0.1),spec,x=x_cif)
spec$t_cif <- t_cif
spec$cif_xs <- cif_xs
spec$x_cif <- x_cif

if (!dir.exists(paste0(folder,"spec"))) dir.create(paste0(folder,"spec"), recursive = TRUE)
saveRDS(spec, paste0(folder,"spec/spec",specnum,".rds"))
