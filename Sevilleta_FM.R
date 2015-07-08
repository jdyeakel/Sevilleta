library(Rcpp)
library(plotrix)
library(RColorBrewer)
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
max_enc <- 5
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
  a <- 1,
  b <- -0.25,
  theta_max <- 50,
  tmax <- 50,
  pk <- pk,
  gain <- gain
)

W <- Cout[[1]]
jstar <- Cout[[2]]
dec <- Cout[[3]]

pal <- brewer.pal(5,"Set1")
resources = c("NA","C3 veg", "C3 seeds", "C4 veg", "C4 seeds", "Insects")
theta <- 10
xx <- jstar[[theta]] +1
pal.m <- as.character(xx); 
pal.m[which(pal.m == "1")] = "white"; pal.m[which(pal.m == "2")] = pal[1]; pal.m[which(pal.m == "3")] = pal[2]; 
pal.m[which(pal.m == "4")] = pal[3]; pal.m[which(pal.m == "5")] = pal[4]; pal.m[which(pal.m == "6")] = pal[5]
lbs <- unique(as.numeric(xx))
par(mar=c(5,5,2,10))
color2D.matplot(xx,extremes=lbs, border=NA, axes=TRUE, xlab="Time", ylab="Energetic Reserves",main=paste("Theta = ",theta-1),cellcolors = pal.m)
legend(tmax,Mc*10,legend=resources[sort(lbs)],pch=22,pt.bg=c("white",pal)[sort(lbs)],xpd=TRUE, bty="n")

