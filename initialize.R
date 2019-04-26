source('readdata.R')

StartYear <- Kayadata[1,1]
NHistYear <- nrow(Kayadata)
EndHYear <- Kayadata[NHistYear,1]
NYear <- NHistYear+82
CO2emit <- matrix(NA,NYear,1)
population <- matrix(NA,NYear,1)
dpop <- matrix(0,NYear,1)
ddpop <- matrix(0,NYear,1)
income <- matrix(NA,NYear,1)
energyint <- matrix(NA,NYear,1)
carbonint <- matrix(NA,NYear,1)
population[1:NHistYear] <- Kayadata[1:NHistYear,2]
pophist <- population
income[1:NHistYear] <- Kayadata[1:NHistYear,3]
inchist <- income
energyint[1:NHistYear] <- Kayadata[1:NHistYear,4]
eihist <- energyint
carbonint[1:NHistYear] <- Kayadata[1:NHistYear,5]
cihist <- carbonint
for (t in 1:NHistYear){
  CO2emit[t] <- population[t]*income[t]*energyint[t]*carbonint[t]/1000000000
}
co2hist <- CO2emit

for (t in 2:NHistYear){
  dpop[t] <- population[t]/population[t-1] -1
}

for (t in 3:NHistYear){
  ddpop[t] = dpop[t] - dpop[t-1]
}

addpop <- 0
for (t in (NHistYear-50):NHistYear){
  addpop <- addpop + ddpop[t]
}
addpop <- addpop/50

Year <- matrix(0,NYear,1)
Year[1:NHistYear] <- Kayadata[1:NHistYear,1]

Temperature <- matrix(NA,NYear,1)
Temperature[1:NHistYear] <- Tempdata[1:NHistYear,2]

Impact <- matrix(NA,NYear,1)
ImpactP <- matrix(NA,NYear,1)

source('initcarboncycle.R')
source('initclimate.R')