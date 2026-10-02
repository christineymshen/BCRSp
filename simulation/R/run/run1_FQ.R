# catching errors

specnum <- 1
runnum <- 1
library(survival)

runlabel <- "FQ"
nsim <- 1000

folder <- "simulation/"
source(paste0(folder,"R/functions.R"))
outfolder <- paste0(folder,"res/", runlabel)
if (!file.exists(outfolder)) dir.create(outfolder, recursive = T)

spec <- readRDS(paste0(folder, "spec/spec",specnum,".rds"))

set.seed(1)

data_kmeans <- kmeans(spec$d,centers = spec$ncenters)

if (spec$isModels){
  X <- cbind(spec$X,comorbidity=spec$W)
  p <- spec$p+1
} else {
  X <- spec$X
  p <- spec$p
}

label1 <- paste(colnames(X),collapse = "+")

df_group <- data.frame(group=data_kmeans$cluster) %>% 
  mutate(group=factor(group))

X2 <- model.matrix(~group-1, data=df_group)[,-1]
X2 <- apply(X2,2,scale,scale=F)
X3 <- cbind(X,X2)
label2 <- paste(colnames(X3),collapse = "+")

i <- j <- 1
idx <- numeric(nsim)

while (i<=nsim){
  
  output_name <- paste0(outfolder,"/spec",specnum,"_run",runnum,"_",i,".rds")
  
  sim <- simdata10(spec,seed=j)
  
  data1 <- cbind(sim$data,X) %>%
    mutate(event=factor(event))
  
  data2 <- cbind(data1,X2)
  
  
  mod1 <- tryCatch({
    coxph(as.formula(paste0("Surv(time,event) ~",label1)), data1, id=c(1:spec$n),
          control=coxph.control(iter.max=100))
  },
  error = function(e) {
    return(NULL)
  })
  
  if (!is.null(mod1)){
    mod2 <- tryCatch({
      coxph(as.formula(paste0("Surv(time,event) ~",label2)), data2, id=c(1:spec$n),
            control=coxph.control(iter.max=100))
    },
    error = function(e) {
      return(NULL)
    })
  } else {mod2 <- NULL}
    
  if (!is.null(mod1) & !is.null(mod2)){
    output <- list(coxph1_m=mod1$coefficients,
                   coxph1_sd=sqrt(diag(mod1$var)),
                   coxph2_m=mod2$coefficients[c(1:p,(p+spec$ncenters):(2*p+spec$ncenters-1))],
                   coxph2_sd=sqrt(diag(mod2$var))[c(1:p,(p+spec$ncenters):(2*p+spec$ncenters-1))])
    
    saveRDS(output, output_name)
    idx[i] <- j
    i <- i + 1
  }
  
  j <- j+1
    
}

saveRDS(idx, paste0(outfolder,"/spec",specnum,"_run",runnum,"_seed.rds"))

