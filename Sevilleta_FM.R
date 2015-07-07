library(Rcpp)
sourceCpp("src/SDP_beq_func.cpp")


#Sevilleta foraging model
#Setup
num_res <- 5

# 1) C3 Veg
# 2) C3 Seeds
# 3) C4 Veg
# 4) C4 Seeds
# 5) Insects

#Probability of finding k resources of food type j
m <- c(4,3,5,2,1)
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
gain <- c(2,3,1,3,5)

#Run SDP
Cout <- SDP_beq_func(
  Mc <- 5,
  a <- 10,
  b <- -0.25,
  theta_max <- 10,
  tmax <- 10,
  pk <- pk,
  gain <- gain
)
