#Analysis of rainfall data


year <- seq(1,365)

monstart <- c(
  181,
  183,
  195,
  195,
  178,
  176,
  168,
  178,
  177,
  207,
  185,
  203
  )

hist(monstart,xlim=c(0,365),breaks=5,col="gray")