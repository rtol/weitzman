ConvF <- 0.00047
BoxShare = matrix(c(0.13, 0.20, 0.32, 0.25, 0.10),nrow=1,ncol=5)
BoxLife = matrix(c(0, 1-exp(-1/363), 1-exp(-1/74), 1-exp(-1/17), 1-exp(-1/2)),nrow=1,ncol=5)

for (t in 2:NHistYear){
  for (b in 1:NCO2Box){
    CO2box[t,b] = CO2box[t-1,b]*(1-BoxLife[b]) + BoxShare[b]*ConvF*CO2emit[t-1,2]
  }
}