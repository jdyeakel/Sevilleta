rm(list=c(ls()))


library(Rcpp)
library(plotrix)
library(RColorBrewer)
#sourceCpp("src/SDP_beq_func.cpp")
sourceCpp("src/SDP_beq_func_trim.cpp")
source("R/timetotheta.R")

colors <- brewer.pal(5,"Spectral")

#Sevilleta foraging model
#Setup
num_res <- 5

# 1) C3 Veg
# 2) C3 Seeds
# 3) C4 Veg
# 4) C4 Seeds
# 5) Insects



#Probability of finding k resources of food type j
m <- c(6.5,0.8,0.4,0.05,2.2) #Good Winter
#m <- c(1.5,0.2,5.3,0.6,2.2) #Good Monsoon
#nu <- c(7,4,10,5,2)
nu <- c(1,1,1,1,1)
max_enc <- 10
pk <- matrix(0,(max_enc+1),num_res)
for (j in 1:num_res) {
  for (k in 1:(max_enc+1)) {
    pk[k,j] <-  dnbinom(k,mu = m[j], size = nu[j])
  }
  pk[,j] <- pk[,j] / sum(pk[,j])
}
#plot negative binomial distributions
plot(seq(0,max_enc,1),pk[,1],type="l",col=colors[1],lwd=3,ylim=c(0,1),
     xlab = "Num. resources encountered", ylab = "Probability")
for (i in 2:num_res) {
  lines(seq(0,max_enc,1),pk[,i],col=colors[i],lwd=3)
}
resources = c("C3 veg", "C3 seeds", "C4 veg", "C4 seeds", "Insects")
legend(max_enc-4,1,legend=resources,pch=22,pt.bg=colors,xpd=TRUE, bty="n")
meank <- numeric(num_res)
for (i in 1:num_res) {
  meank[i] <- seq(0,max_enc) %*% pk[,i]
}

#Resource gain
#each unit is 10 kJ/gram
gain <- c(1.5,2.1,1.5,2.1,2.5)

#Run SDP
Cout <- SDP_beq_func_trim(
  Mc <- 40,
  a <- 0.2,
  b <- 0.75,
  theta_max <- 10*Mc,
  tmax <- 100,
  pk <- pk,
  gain <- gain
)

W <- Cout[[1]]
jstar <- Cout[[2]]
dec <- Cout[[3]]
ttt <- timetotheta(W,jstar)
W_xt <- ttt[[1]]; jstar_xt <- ttt[[2]]


###########
# PLOTTING
###########
#Find the 'eat cache' options
xxjstar <- jstar
cache_threshold <- 0.8
for (t in 1:(tmax-1)) {
  for (theta in 1:theta_max) {
    for (x in (floor(0.5*Mc)+1):Mc) {
      j = jstar[[t]][x,theta]
      #dec_vector = dec[[t]][[theta]][[x]][,j]
      prob_none = pk[1,j]
      if (prob_none > cache_threshold) {
        xxjstar[[t]][x,theta] <- 6
      }
    }
  }
}
#pal <- brewer.pal(5,"Set1")
par(mfrow=c(1,3))
timeseq <- c(1,75,95)
tic <- 0
for (t in timeseq) {
  tic <- tic + 1
  pal <- colors
  resources = c("NA","C3 veg", "C3 seeds", "C4 veg", "C4 seeds", "Insects","Pr(>80%) Cache")
  time <- t
  xx <- xxjstar[[time]] +1
  pal.m <- as.character(xx); 
  pal.m[which(pal.m == "1")] = "white"; pal.m[which(pal.m == "2")] = pal[1]; pal.m[which(pal.m == "3")] = pal[2]; 
  pal.m[which(pal.m == "4")] = pal[3]; pal.m[which(pal.m == "5")] = pal[4]; pal.m[which(pal.m == "6")] = pal[5]
  pal.m[which(pal.m == "7")] = "gray"
  lbs <- unique(as.numeric(xx))
  #par(mar=c(5,5,2,10))
  color2D.matplot(xx,extremes=lbs, border=NA, axes=TRUE, xlab="Cache reserves (1 unit = 10 kJ)", ylab="Energetic reserves (1 unit = 10 kJ",main=paste("Time = ",time),cellcolors = pal.m)
  if (tic == 1) {
    legend(0,Mc,legend=resources[sort(lbs)],pch=22,pt.bg=c("white",pal,"gray")[sort(lbs)],xpd=TRUE, bty="n")
  }
}


## TTT Plots
pal <- colors
resources = c("NA","C3 veg", "C3 seeds", "C4 veg", "C4 seeds", "Insects")
theta <- 100
xx <- jstar_xt[[theta]] +1
pal.m <- as.character(xx); 
pal.m[which(pal.m == "1")] = "white"; pal.m[which(pal.m == "2")] = pal[1]; pal.m[which(pal.m == "3")] = pal[2]; 
pal.m[which(pal.m == "4")] = pal[3]; pal.m[which(pal.m == "5")] = pal[4]; pal.m[which(pal.m == "6")] = pal[5]
lbs <- unique(as.numeric(xx))
par(mar=c(5,5,2,10))
color2D.matplot(xx,extremes=lbs, border=NA, axes=TRUE, xlab="Time", ylab="Energetic Reserves",main=paste("Theta = ",theta-1),cellcolors = pal.m)
legend(tmax,Mc,legend=resources[sort(lbs)],pch=22,pt.bg=c("white",pal,"gray")[sort(lbs)],xpd=TRUE, bty="n")








