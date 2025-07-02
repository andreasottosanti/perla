library(perla)
set.seed(123)

N <- 20
x <- list()
for(i in 1:N){
  print(i)
  x[[i]] <- generate.simulations(spatial.map = west_states_data,
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


for(i in 1:20){
  res <- perla(y = x[[i]]$map, K = 4, rho.value = 0.99, mean.penalty = c("c"))
}
