rm(list = ls())
# Running models ----------------------------------------------------------
R <- 10^2
burnin <- 1:50
N <- 1#length(dir("~/perla/Simulations/Simulation1/Data/"))

# Penalisation c ----------------------------------------------------------
Rand_values <- DIC_values <- numeric(N)
for(i in 1:N){
  values <- list()
  x <- readRDS(paste("~/perla/Simulations/Simulation1/Data/Data",i,".RDS",sep=""))
  for(k in 1:4)
    values[[k]] <- perla(y = x$map, K = 4, R = R, burnin = burnin,
                         mean.penalty = c("c"),
                         rho.value = 0.99,
                         initialization = list(Sigma = x$Sigma), seed = 123*k)
  results <- merge.perla(values)
  results <- recover.loglikelihood(value = results)
  relab <- remove.label.switching(values = results, loglikelihood.values = results$loglik)
  results <- perla::information.criteria(values=results)
  results$rand <- mclust::adjustedRandIndex(x$z, relab$relabelling$clusters[2,])
  save(results, relab, file = paste("~/perla/Simulations/Simulation1/Results/c/Results",i,".RData",sep=""))
}

# Penalisation d ----------------------------------------------------------
Rand_values <- DIC_values <- numeric(N)
for(i in 1:N){
  values <- list()
  x <- readRDS(paste("~/perla/Simulations/Simulation1/Data/Data",i,".RDS",sep=""))
  for(k in 1:4)
    values[[k]] <- perla(y = x$map, K = 4, R = R, burnin = burnin,
                         mean.penalty = c("d"),
                         rho.value = 0.99,
                         initialization = list(Sigma = x$Sigma), seed = 123*k)
  results <- merge.perla(values)
  results <- recover.loglikelihood(value = results)
  relab <- remove.label.switching(values = results, loglikelihood.values = results$loglik)
  results <- perla::information.criteria(values=results)
  results$rand <- mclust::adjustedRandIndex(x$z, relab$relabelling$clusters[2,])
  save(results, relab, file = paste("~/perla/Simulations/Simulation1/Results/d/Results",i,".RData",sep=""))
}

# Penalisation c,d ----------------------------------------------------------
Rand_values <- DIC_values <- numeric(N)
for(i in 1:N){
  values <- list()
  x <- readRDS(paste("~/perla/Simulations/Simulation1/Data/Data",i,".RDS",sep=""))
  for(k in 1:4)
    values[[k]] <- perla(y = x$map, K = 4, R = R, burnin = burnin,
                         mean.penalty = c("c","d"),
                         rho.value = 0.99,
                         initialization = list(Sigma = x$Sigma), seed = 123*k)
  results <- merge.perla(values)
  results <- recover.loglikelihood(value = results)
  relab <- remove.label.switching(values = results, loglikelihood.values = results$loglik)
  results <- perla::information.criteria(values=results)
  results$rand <- mclust::adjustedRandIndex(x$z, relab$relabelling$clusters[2,])
  save(results, relab, file = paste("~/perla/Simulations/Simulation1/Results/c_d/Results",i,".RData",sep=""))
}

# Penalisation cd ----------------------------------------------------------
Rand_values <- DIC_values <- numeric(N)
for(i in 1:N){
  values <- list()
  x <- readRDS(paste("~/perla/Simulations/Simulation1/Data/Data",i,".RDS",sep=""))
  for(k in 1:4)
    values[[k]] <- perla(y = x$map, K = 4, R = R, burnin = burnin,
                         mean.penalty = c("cd"),
                         rho.value = 0.99,
                         initialization = list(Sigma = x$Sigma), seed = 123*k)
  results <- merge.perla(values)
  results <- recover.loglikelihood(value = results)
  relab <- remove.label.switching(values = results, loglikelihood.values = results$loglik)
  results <- perla::information.criteria(values=results)
  results$rand <- mclust::adjustedRandIndex(x$z, relab$relabelling$clusters[2,])
  save(results, relab, file = paste("~/perla/Simulations/Simulation1/Results/cd/Results",i,".RData",sep=""))
}

# Penalisation d,cd ----------------------------------------------------------
Rand_values <- DIC_values <- numeric(N)
for(i in 1:N){
  values <- list()
  x <- readRDS(paste("~/perla/Simulations/Simulation1/Data/Data",i,".RDS",sep=""))
  for(k in 1:4)
    values[[k]] <- perla(y = x$map, K = 4, R = R, burnin = burnin,
                         mean.penalty = c("d","cd"),
                         rho.value = 0.99,
                         initialization = list(Sigma = x$Sigma), seed = 123*k)
  results <- merge.perla(values)
  results <- recover.loglikelihood(value = results)
  relab <- remove.label.switching(values = results, loglikelihood.values = results$loglik)
  results <- perla::information.criteria(values=results)
  results$rand <- mclust::adjustedRandIndex(x$z, relab$relabelling$clusters[2,])
  save(results, relab, file = paste("~/perla/Simulations/Simulation1/Results/d_cd/Results",i,".RData",sep=""))
}

# Penalisation NULL ----------------------------------------------------------
Rand_values <- DIC_values <- numeric(N)
for(i in 1:N){
  values <- list()
  x <- readRDS(paste("~/perla/Simulations/Simulation1/Data/Data",i,".RDS",sep=""))
  for(k in 1:4)
    values[[k]] <- perla(y = x$map, K = 4, R = R, burnin = burnin,
                         mean.penalty = NULL,
                         rho.value = 0.99,
                         initialization = list(Sigma = x$Sigma), seed = 123*k)
  results <- merge.perla(values)
  results <- recover.loglikelihood(value = results)
  relab <- remove.label.switching(values = results, loglikelihood.values = results$loglik)
  results <- perla::information.criteria(values=results)
  results$rand <- mclust::adjustedRandIndex(x$z, relab$relabelling$clusters[2,])
  save(results, relab, file = paste("~/perla/Simulations/Simulation1/Results/NI/Results",i,".RData",sep=""))
}




