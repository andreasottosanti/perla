rm(list = ls())
library(perla)

saving_directory <- "/mnt/callisto/Sottosanti/perla/Applications/US_East_Coast/Results"

# Running models ----------------------------------------------------------
R <- 10^4
burnin <- 1:5000
gender <- c("males", "females")
save_penalisation <- c("c", "d", "c_d", "cd", "d_cd")
code_penalisation <- list(c("c"),
                          c("d"),
                          c("c","d"),
                          c("cd"),
                          c("d","cd"))

for(i in gender){                                                                     # run on males and females
  x <- readRDS(paste("~/perla/Applications/US_East_Coast/Data/",i,".RDS",sep=""))
  W <- readRDS("~/perla/Applications/US_East_Coast/Data/proximity_matrix.RDS")
  x@data <- log(x@data[,6:8])
  if(i == "males") K_to_try <- 2:3
  if(i == "females") K_to_try <- 2:4
  for(p in 1:length(save_penalisation)){
    for(k in K_to_try){
      values <- list()
      for(j in 1:4){
        values[[j]] <- perla(y = x, W = W,
                             K = k, R = R, burnin = burnin,
                             mean.penalty = code_penalisation[[p]],
                             rho.value = 0.99,
                             initialization = list(Sigma = cov(x@data)), seed = 1234*j)
      }
      results <- merge.perla(values)
      results <- recover.loglikelihood(value = results)
      relab <- remove.label.switching(values = results, loglikelihood.values = results$loglik)
      results <- perla::information.criteria(values=results)
      if(!(i %in% dir(saving_directory))) dir.create(paste(saving_directory,i,sep="/"))
      save(results, relab,
           file = paste(saving_directory,"/",i,"/Results_",k,"_",save_penalisation[p],".RData",sep=""))
    }
  }
}



# females
# 124
i <- "females"
x <- readRDS(paste("~/perla/Applications/US_East_Coast/Data/",i,".RDS",sep=""))
W <- readRDS("~/perla/Applications/US_East_Coast/Data/proximity_matrix.RDS")
x@data <- log(x@data[,6:8])
if(i == "males") K_to_try <- 2:3
if(i == "females") K_to_try <- 2:4
for(p in 1:length(save_penalisation)){
  for(k in K_to_try){
    values <- list()
    for(j in 1:4){
      values[[j]] <- perla(y = x, W = W,
                           K = k, R = R, burnin = burnin,
                           mean.penalty = code_penalisation[[p]],
                           rho.value = 0.99,
                           initialization = list(Sigma = cov(x@data)), seed = 129*j)
    }
    results <- merge.perla(values)
    results <- recover.loglikelihood(value = results)
    relab <- remove.label.switching2(values = results, loglikelihood.values = results$loglik)
    results <- perla::information.criteria(values=results)
    if(!(i %in% dir(saving_directory))) dir.create(paste(saving_directory,i,sep="/"))
    save(results, relab,
         file = paste(saving_directory,"/",i,"/Results_",k,"_",save_penalisation[p],".RData",sep=""))
  }
}
