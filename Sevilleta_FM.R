library(Rcpp)
library(plotrix)
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
max_enc <- 30
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
  a <- 4,
  b <- -0.25,
  theta_max <- 100,
  tmax <- 200,
  pk <- pk,
  gain <- gain
)

W <- Cout[[1]]
jstar <- Cout[[2]]
dec <- Cout[[3]]

xx <- jstar[[100]]
lbs <- unique(as.numeric(xx))
par(mar=c(3,3,1,10))
color2D.matplot(xx,extremes=lbs+1, border=NA, axes=TRUE, xlab="", ylab="",main="")
legend(tmax,50,legend=as.character(lbs+1),pch=22,pt.bg=lbs+1,xpd=TRUE, bty="n")

