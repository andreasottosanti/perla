library(perla)
set.seed(123)


# Data generation ---------------------------------------------------------
N <- 20
for(i in 1:N){
  print(i)
  x <- generate.simulations(spatial.map = west_states_data,
                            K = 3,
                            d = 10,
                            Sigma = NULL,
                            range.mu = 0.5,
                            range.Sigma = 0.1,
                            rho = runif(2, 0.8, 1),
                            prob.null.centroid = 0.5,
                            scale.factor.variance = 0.05,
                            correct.mean.prob = F,
                            plot.map = F)
  saveRDS(x, file = paste("Simulations/Simulation1/Data/Data",i,".RDS", sep = ""))
}



# Running models ----------------------------------------------------------
Rand_values <- DIC_values <- numeric(N)
for(i in 1:N){
  values <- list()
  x <- readRDS(paste("~/perla/Simulations/Simulation1/Data/Data",i,".RDS",sep=""))
  for(k in 1:4)
    values[[k]] <- perla(y = x$map, K = 4, R = 10^2, burnin = 1:50,
                         mean.penalty = c("c"),
                         rho.value = 0.99,
                         initialization = list(Sigma = x$Sigma), seed = 123*k)
  results <- merge.perla(values)
  results <- recover.loglikelihood(value = results)
  relab <- remove.label.switching(values = results, loglikelihood.values = results)
  results <- perla::information.criteria(values=results)
  results$rand <- mclust::adjustedRandIndex(x$z, relab$relabelling$clusters[2,])
  save(results, rand, DIC, loglikelihood.values, relab, file = paste("~/perla/Simulations/Simulation1/Results/Results",i,".RDS",sep=""))
}
save(results, rand, DIC, loglikelihood.values, relab, file = paste("SIM1_results1_c.Rdata",sep=""))
rm(list=ls())

for(i in 1:20){
  res <- perla(y = x[[i]]$map, K = 4, rho.value = 0.99, mean.penalty = c("c"))
}
