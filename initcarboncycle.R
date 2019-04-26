NCO2Box <- 5
CO2box <- matrix(0,NYear,NCO2Box)
CO2pre <- 275
CO2box[1,1] <- CO2pre
CO2boxP <- CO2box
CO2conc <- matrix(0,NYear,1)
CO2conc[1] <- CO2pre
CO2concP <- CO2conc
ConvF <- 0.000471
BoxShare = matrix(c(0.13, 0.20, 0.32, 0.25, 0.10),nrow=1,ncol=5)
BoxLife = matrix(c(0, 1-exp(-1/363), 1-exp(-1/74), 1-exp(-1/17), 1-exp(-1/2)),nrow=1,ncol=5)