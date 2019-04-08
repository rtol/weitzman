source('readdata.R')
StartYear <- CO2data[1,1]
NHistYear <- nrow(CO2data)
EndHYear <- CO2data[NHistYear,1]
NYear <- NHistYear+82
CO2emit <- matrix(0,NYear,1)
CO2emit[1:NHistYear] <- CO2data[1:NHistYear,2]
Year <- matrix(0,NYear,1)
Year[1:NHistYear] <- CO2data[1:NHistYear,1]

Temperature <- matrix(0,NYear,1)
Temperature[1:NHistYear] <- Tempdata[1:NHistYear,2]

source('initcarboncycle.R')
source('initclimate.R')