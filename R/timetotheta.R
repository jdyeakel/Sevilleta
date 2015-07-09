timetotheta <- function(W,jstar) {
  tmax <- length(W)
  thetamax <- length(W[[1]][1,])
  xmax <- length(W[[1]][,1])
  
  jstar_xt <- list()
  W_xt <- list()
  for (theta in 1:theta_max) {
    jstarm <- matrix(0,xmax,tmax-1)
    Wm <- matrix(0,xmax,tmax-1)
    for (x in 1:xmax) {
      for (t in 1:(tmax-1)) {
        jstarm[x,t] <- jstar[[t]][x,theta]
        Wm[x,t] <- W[[t]][x,theta]
      }
    }
    jstar_xt[[theta]] <- jstarm
    W_xt[[theta]] <- Wm
  }
  
  Rout <- list()
  Rout[[1]] <- W_xt
  Rout[[2]] <- jstar_xt

  return(Rout) 
  
}