library(Rcpp)
sourceCpp("src/SDP_beq_func.cpp")


#Sevilleta foraging model
#Setup
num_res <- 5

#Probability of finding k resources of food type j
m <- c(1,2,3,4,5)
nu <- c(1,1,1,1,1)
max_enc <- 50
pk <- matrix(0,(max_enc+1),num_res)
for (j in 1:num_res) {
  for (k in 1:(max_enc+1)) {
    pk[k,j] <-  dnbinom(k,mu = m[j], size = nu[j])
  }
  pk[,j] <- pk[,j] / sum(pk[,j])
}

#Resource gain
gain <- c(1,2,3,4,5)
#Resource cost
cost <- c(1,2,3,4,5)

#Run SDP
Cout <- SDP_beq_func(
  Mc < -5,
  a <- 0.1,
  b <- 0.75,
  tmax <- 10,
  pk <- pk
  gain <- gain,
  cost <- cost
)
