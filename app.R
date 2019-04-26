# Shiny Integrated Assessment Model

library(shiny)

# Define UI for user input
ui <- fluidPage(
  titlePanel("Hit the Paris target: Keep temperature below 2 degrees! (Don't move too far from the data.)"),
  fluidRow(
    column(3,
      sliderInput("ClimSens",
        "Warming per CO2 doubling",
        min = 0.1,
        max = 10,
        value = 3.9),
      helpText("Equilibrium warming for doubling of atmospheric CO2."),
      sliderInput("PopRate",
        "Rate population growth falls",
        min = -0.118,
        max = 0.079,
        value = -0.019),
      helpText("The decline, in percentage points, of the population growth rate."),
      sliderInput("IncRate",
        "Income growth",
        min = -2.6,
        max = 5.8,
        value = 1.6),
      helpText("The growth rate of per capita income, in percentage."),
      sliderInput("AEEI",
        "Energy efficiency improvement",
        min = -2.0,
        max = 3.6,
        value = 0.8),
      helpText("The rate of energy efficiency of improvement, in percentage."),
      sliderInput("ACEI",
        "Carbon efficiency improvement",
        min = -2.9,
        max = 3.3,
        value = 0.2),
      helpText("The rate of decarbonisation of the energy supply, in percentage."),
      sliderInput("Topt",
        "Optimal temperature",
        min = -1.0,
        max = 5.0,
        step = 0.01,
        value = -0.38),
      helpText("The temperature that is best for humans, difference to pre-industrial."),
      sliderInput("Benchmark",
        "Impact of 2.5K warming",
        min = -10.0,
        max = 5.0,
        step = 0.01,
        value = -1.34),
      helpText("The welfare impact of 2.5K warming, in percent of income."),
      sliderInput("IncElas",
        "Income elasticity of impact",
        min = -3.5,
        max = 1.0,
        step = 0.01,
        value = -1.68),
      helpText("The change in vulnerability due to a change in income."),
      sliderInput("DiscRate",
        "Discount rate",
        min = 0.0,
        max = 10.0,
        step = 0.1,
        value = 5.0),
      helpText("The consumption discount rate, in percent per year."),
      tags$a(href="https://github.com/rtol/SIAM", "Click here for code and data on GitHub")
    ),
    column(3,
      plotOutput("plotP"),
      plotOutput("plotCI"),
      plotOutput("plotT"),
      textOutput("SCC")
    ),
    column(3,
      plotOutput("plotY"),
      plotOutput("plotE"),
      plotOutput("plotI"),
      helpText("All results are either global totals or global averages. circles = observations, lines = projections.")
      ),
    column(3,
      plotOutput("plotEI"),
      plotOutput("plotC"),
      plotOutput("plotD"),
      helpText("Likelihood ratio is relative to the best fit to the 1969-2018 observations (1850-2018 for temperature).")
    )
  )
)

# Define server to generate output
server <- function(input, output) {
  
  resultsout <- reactive({
    
    source("initialize.R")
    
    RFpar <- input$ClimSens/l2/log(2)
    RFpar0 <- 3.9/l2/log(2)
    popgr <- dpop[NHistYear]
    addpop <- input$PopRate/100
    ppop <- matrix(exp(-0.5*((addpop+0.00019)/0.00033)^2),NYear,1)
    
    incgr <- input$IncRate/100
    pinc <- matrix(exp(-0.5*((incgr-0.016)/0.014)^2),NYear,1)
    
    eigr <- -input$AEEI/100
    pei <- matrix(exp(-0.5*((eigr+0.008)/0.009)^2),NYear,1)
    
    cigr <- -input$ACEI/100
    pci <- matrix(exp(-0.5*((cigr+0.002)/0.010)^2),NYear,1)
    
    beta <- input$Benchmark/(2.5*2.5-2*2.5*input$Topt)
    alpha <- -input$Topt*2*beta
    
    SCC <- 0
    
    for (t in 2:NYear){
      if (t > NHistYear){
        Year[t] <- Year[t-1]+1
        popgr <- popgr + addpop
        population[t] <- population[t-1]*(1+popgr)
        income[t] <- income[t-1]*(1+incgr)
        energyint[t] <- energyint[t-1]*(1+eigr)
        carbonint[t] <- carbonint[t-1]*(1+cigr)
        CO2emit[t] <- population[t]*income[t]*energyint[t]*carbonint[t]/1000000000
      }
      for (b in 1:NCO2Box){
        CO2box[t,b] <- CO2box[t-1,b]*(1-BoxLife[b]) + BoxShare[b]*ConvF*CO2emit[t-1]
      }
      for (b in 1:NCO2Box){
        if (t == NHistYear+2){
          CO2boxP[t,b] <- CO2boxP[t-1,b]*(1-BoxLife[b]) + BoxShare[b]*ConvF*(CO2emit[t-1]+1)
        }
        else {
          CO2boxP[t,b] <- CO2boxP[t-1,b]*(1-BoxLife[b]) + BoxShare[b]*ConvF*CO2emit[t-1]
        }
      }
      CO2conc[t] <- sum(CO2box[t,1:NCO2Box])
      RadForc[t] <- RFpar*log(CO2conc[t]/CO2pre)
      TempAtm[t] <- TempAtm[t-1] + l1*(l2*RadForc[t]-TempAtm[t-1]) + l3*(TempOc[t-1]-TempAtm[t-1])
      TempOc[t] <- TempOc[t-1] + l4*(TempAtm[t-1]-TempOc[t-1])
      RadForc0[t] <- RFpar0*log(CO2conc[t]/CO2pre)
      TempAtm0[t] <- TempAtm0[t-1] + l1*(l2*RadForc0[t]-TempAtm0[t-1]) + l3*(TempOc0[t-1]-TempAtm0[t-1])
      TempOc0[t] <- TempOc0[t-1] + l4*(TempAtm0[t-1]-TempOc0[t-1])
      Impact[t] <- (alpha*TempAtm[t] + beta*TempAtm[t]^2)*(income[t]/income[NHistYear])^input$IncElas
      CO2concP[t] <- sum(CO2boxP[t,1:NCO2Box])
      RadForcP[t] <- RFpar*log(CO2concP[t]/CO2pre)
      TempAtmP[t] <- TempAtmP[t-1] + l1*(l2*RadForcP[t]-TempAtmP[t-1]) + l3*(TempOcP[t-1]-TempAtmP[t-1])
      TempOcP[t] <- TempOcP[t-1] + l4*(TempAtmP[t-1]-TempOcP[t-1])
      ImpactP[t] <- (alpha*TempAtmP[t] + beta*TempAtmP[t]^2)*(income[t]/income[NHistYear])^input$IncElas
      if (t <= NHistYear){
        DF <- 0
      }
      else if (t == NHistYear+1){
        DF <- 1
      }
      else {
        DF <- DF/(1+input$DiscRate/100)
      }
      SCC <- SCC + (Impact[t] - ImpactP[t])*income[t]*population[t]*DF
    }
    
    #now run the model backwards in time
    popgr <- dpop[NHistYear]
    for (t in (NHistYear-1):1){
      popgr <- max(0,popgr - addpop)
      population[t] <- population[t+1]/(1+popgr)
      income[t] <- income[t+1]/(1+incgr)
      energyint[t] <- energyint[t+1]/(1+eigr)
      carbonint[t] <- carbonint[t+1]/(1+cigr)
      co2hist[t] <- population[t]*income[t]*energyint[t]*carbonint[t]/1000000000
    }
    
    resid <- TempAtm - Temperature
    resid0 <- TempAtm0 - Temperature
    RSS <- 0
    RSS0 <- 0
    mtemp <- 0
    mtsq <- 0
    n <- 0
    for (t in 101:NHistYear){
      n <- n+1
      RSS <- RSS + (resid[t]-0.3)^2
      RSS0 <- RSS0 + (resid0[t]-0.3)^2
      mtemp <- mtemp + Temperature[t]-0.3
      mtsq <- mtsq + (Temperature[t]-0.3)^2
    }
    mtemp <- mtemp/n
    mtsq <- mtsq/n
    vartemp <- mtsq-mtemp*mtemp
    p <- exp(-0.5*RSS/vartemp)
    p0 <- exp(-0.5*RSS0/vartemp)
    ptmp <- matrix(p/p0,NYear,1)
    SCCv <- matrix(SCC,NYear,1)
    #                     1     2          3       4         5           6          7       8          9          10  11     12   13  14     15        16         17     18     19        20    21                
    resultsout <- cbind(Year, CO2emit, CO2conc, TempAtm, Temperature, population, income, energyint, carbonint, ppop, pinc, pei, pci, ptmp, pophist, inchist, eihist, cihist, co2hist, ImpactP, SCCv)
  })
  
    resultsimp <- reactive({
      impdata <- read.csv(file="totalimpact.csv", header=FALSE, sep=",")
      temp <- impdata[1:27,1]
      impobs <- impdata[1:27,2]
      beta <- input$Benchmark/(2.5*2.5-2*2.5*input$Topt)
      alpha <- -input$Topt*2*beta
      impmod <- matrix(NA,27,1)
      impmod1 <- matrix(NA,27,1)
      impbest <- matrix(NA,27,1)
      res <- matrix(NA,27,1)
      res0 <- matrix(NA,27,1)
      tempgrid <- matrix(NA,27,1)
      tempgrid[1] <- -1
      for (i in 2:27){
        tempgrid[i] <- tempgrid[i-1] + 7/26
      }
      
      for (i in 1:27){
        impmod[i] <- alpha*tempgrid[i] + beta*tempgrid[i]^2
        impmod1[i] <- alpha*temp[i] + beta*temp[i]^2
        impbest[i] <- -0.12496*temp[i] -0.16442*temp[i]^2
        res[i] <- impobs[i] - impmod1[i]
        res0[i] <- impobs[i] - impbest[i]
      }
      
      rss <- 0
      rss0 <- 0
      mimp <- 0
      misq <- 0
      for (i in 1:27){
        rss <- rss + res[i]^2
        rss0 <- rss0 + res0[i]^2
        mimp <- mimp + impobs[i]
        misq <- misq + impobs[i]^2
      }
      mimp <- mimp/27
      misq <- misq/27
      varimp <- misq-mimp*mimp
      pi <- exp(-0.5*rss/varimp)
      pi0 <- exp(-0.5*rss0/varimp)
      pimp <- matrix(pi/pi0,27,1)
      
      resultsimp <- cbind(temp,impobs,tempgrid,impmod,pimp)  
    })
    
    resultsco2 <- reactive({
      resultsco2 <- read.csv(file="co2concentration.csv", header=FALSE, sep=",")    
    })
    
    output$plotI <- renderPlot({
      plot(resultsimp()[1:27,1],resultsimp()[1:27,2], xlab = 'degree Celsius', ylab = 'welfare equivalent income loss')
      lines(resultsimp()[1:27,3],resultsimp()[1:27,4], type="l", lty="solid")
      title(paste("Impact of climate change. Likelihood ratio:",format(resultsimp()[1,5],digits=6)))
    })

    output$plotT <- renderPlot({
      paris <- rep.int(2, NYear)
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,4], type = "l", xlab = 'year', ylab = 'degree Celsius',ylim = c(-0.5,6.0))
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,5]+0.3, type="p")
      lines(resultsout()[1:NYear,1],paris, type="l", lty="dashed")
      title(paste("Temperature. Likelihood ratio:", format(resultsout()[NYear,10]*resultsout()[NYear,11]*resultsout()[NYear,12]*resultsout()[NYear,13]*resultsout()[NYear,14],digits=6)))
    })
    
    output$plotE <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,2]/1000, type = "l", xlab = 'year', ylab = 'billion tonnes of carbon per year')
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,19]/1000, type="p")
      title(paste("Carbon dioxide emissions. Likelihood ratio:",format(resultsout()[NYear,10]*resultsout()[NYear,11]*resultsout()[NYear,12]*resultsout()[NYear,13],digits=6)))
    })
    
    output$plotP <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,6]/1000000000, type = "l", xlab = 'year', ylab = 'billion people')
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,15]/1000000000, type="p")
      title(paste("Population. Likelihood ratio:",format(resultsout()[NYear,10],digits=6)))
    })
    
    output$plotY <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,7], type = "l", xlab = 'year', ylab = 'dollar (2010) per person per year')
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,16], type="p")
      title(paste("Per capita income in 2100. Likelihood ratio:",format(resultsout()[NYear,11],digits=6)))
    })
    
    output$plotEI <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,8]*1000, type = "l", xlab = 'year', ylab = 'gram oil equivalent per dollar', ylim = c(0,300))
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,17]*1000, type="p")
      title(paste("Energy intensity. Likelihood ratio:",format(resultsout()[NYear,12],digits=6)))
    })

    output$plotCI <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,9], type = "l", xlab = 'year', ylab = 'gram carbon per gram oil equivalent')
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,18], type="p")
      title(paste("Carbon intensity. Likelihood ratio:",format(resultsout()[NYear,13],digits=6)))
    })
    
    output$plotC <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,3], type = "l", xlab = 'year', ylab = 'parts per million')
      lines(resultsco2()[,1],resultsco2()[,2], type="p")
      title(paste("Carbon dioxide concentration. Likelihood ratio:",format(resultsout()[NYear,10]*resultsout()[NYear,11]*resultsout()[NYear,12]*resultsout()[NYear,13],digits=6)))
    })
    
    output$plotD <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,20], type = "l", xlab = 'year', ylab = 'percent income')
      title(paste("Impact of climate change. Likelihood ratio:",format(resultsimp()[1,5]*resultsout()[NYear,10]*resultsout()[NYear,11]*resultsout()[NYear,12]*resultsout()[NYear,13]*resultsout()[NYear,14],digits=6)))
    })
    
    output$SCC <- renderText({paste("Social cost of carbon: $",format(resultsout()[NYear,21]/1000000000,digits=4),"/tC")})
}

# Run the application 
shinyApp(ui = ui, server = server)