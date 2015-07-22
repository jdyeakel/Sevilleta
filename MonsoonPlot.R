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

monthstart <- c(1,31,59,90,120,151,181,212,243,273,304,334)
precip <- c(
  7.2,
  6.5,
  15.7,
  11.6,
  8.6,
  10.3,
  37.3,
  30.8,
  22.0,
  22.6,
  10.8,
  12.6)

