# Shiny Integrated Assessment Model

library(shiny)

# Define UI for user input
ui <- fluidPage(
  titlePanel("Compute the social cost of carbon as a function of scenario, impact function, and taste parameters."),
  fluidRow(
    column(3,
      h3(textOutput("SCC")),
      radioButtons("scenario", "Scenario",
                   choices = list("SSP1" = 1, "SSP2" = 2, "SSP3" = 3, "SSP4" = 4, "SSP5" = 5, "Trend" = 6), selected = 6),
      helpText("Scenarios are decomposed according to the Kaya Identity."),
      radioButtons("radio", "Impact model",
           choices = list("Hope (2006)" = 1, "Nordhaus (1991)" = 2, "Nordhaus (2008)" = 3, "Tol (2009)" = 4, "Tol (2019)" = 5, "Weitzman (2012)" = 6),selected = 2),
      helpText("Different authors have proposed different functions to describe the impact of climate change."),
      sliderInput("ClimSens",
        "Warming per CO2 doubling",
        min = 0.1,
        max = 10,
        value = 3.9),
      helpText("Equilibrium warming for doubling of atmospheric CO2."),
      sliderInput("eta",
        "Rate of relative risk aversion",
        min = 0.0,
        max = 10.0,
        step = 0.1,
        value = 1.0),
      helpText("Curvature of the utility function."),
      sliderInput("PRTP",
        "Pure rate of time preference",
        min = 0.0,
        max = 5.0,
        step = 0.1,
        value = 1.0),
      helpText("The utility discount rate, in percent per year."),
      sliderInput("IncElas",
        "Income elasticity of impact",
        min = -3.5,
        max = 1.0,
        step = 0.1,
        value = 0.0),
      helpText("The change in vulnerability due to a change in income."),
      tags$a(href="https://github.com/rtol/SIAM", "Click here for code and data on GitHub")
    ),
    column(3,
      plotOutput("plotT"),
      plotOutput("plotP"),
      plotOutput("plotCI"),
      helpText("All results are global totals or averages.")
    ),
    column(3,
      plotOutput("plotI"),
      plotOutput("plotY"),
      plotOutput("plotE"),
      helpText("circles = observations, lines = projections.")
      ),
    column(3,
      plotOutput("plotD"),
      plotOutput("plotEI"),
      plotOutput("plotC")
    )
  )
)

# Define server to generate output
server <- function(input, output) {
  
  resultsout <- reactive({
    
    source("initialize.R")
    
    RFpar <- input$ClimSens/l2/log(2)
    popgr <- dpop[NHistYear]
    if (input$scenario == "1"){
      addpop <- -0.000285
      incgr <- 0.0239
      eigr <- -0.0192
      cigr <- -0.0058
    }
    else if(input$scenario == "2"){
      addpop <- -0.000208
      incgr <- 0.0203
      eigr <- -0.0131
      cigr <- -0.0022
    }
    else if(input$scenario == "3"){
      addpop <- -0.000110
      incgr <- 0.0089
      eigr <- -0.0053
      cigr <- -0.0003
      
    }
    else if(input$scenario == "4"){
      addpop <- -0.000200
      incgr <- 0.0152
      eigr <- -0.0199
      cigr <- -0.0059
    }
    else if(input$scenario == "5"){
      addpop <- -0.000265
      incgr <- 0.0299
      eigr <- -0.0165
      cigr <- -0.0008
    }
    else{
      addpop <- -0.00011
      incgr <- 0.016
      eigr <- -0.008
      cigr <- -0.002
    }
    
    
    if (input$radio == "1") {
      imp1 <- 0
      imp13 <- -1.0/2.5^1.3
      imp2 <- 0
      imp6 <- 0
    }
    else if (input$radio == "2") {
      imp1 <- 0
      imp13 <- 0
      imp2 <- -1.0/3/3
      imp6 <- 0
    }
    else if (input$radio == "3") {
      imp1 <- 0
      imp13 <- 0
      imp2 <- -2.5/3/3
      imp6 <- 0
    }
    else if (input$radio == "4"){
      imp1 <- 2.46
      imp13 <- 0
      imp2 <- -1.11
      imp6 <- 0
    }
    else if (input$radio == "5"){
      imp1 <- -0.12
      imp13 <- 0
      imp2 <- -0.16
      imp6 <- 0
    }
    else if (input$radio == "6"){
      imp1 <- 0
      imp13 <- 0
      imp2 <- -100/20.46/20.46
      imp6 <- -100/6.081^6.754
      impe <- 0
    }
    else {
      imp1 <- 0
      imp13 <- 0
      imp2 <- 0
      imp6 <- 0
    }
    
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
      Impact[t] <- (imp1*TempAtm[t] + imp13*abs(TempAtm[t])^1.3 + imp2*TempAtm[t]^2 + imp6*abs(TempAtmP[t])^6.754)*(income[t]/income[NHistYear])^input$IncElas
      CO2concP[t] <- sum(CO2boxP[t,1:NCO2Box])
      RadForcP[t] <- RFpar*log(CO2concP[t]/CO2pre)
      TempAtmP[t] <- TempAtmP[t-1] + l1*(l2*RadForcP[t]-TempAtmP[t-1]) + l3*(TempOcP[t-1]-TempAtmP[t-1])
      TempOcP[t] <- TempOcP[t-1] + l4*(TempAtmP[t-1]-TempOcP[t-1])
      ImpactP[t] <- (imp1*TempAtmP[t] + imp13*abs(TempAtmP[t])^1.3 + imp2*TempAtmP[t]^2 + imp6*abs(TempAtmP[t])^6.754)*(income[t]/income[NHistYear])^input$IncElas
      if (t <= NHistYear){
        DF <- 0
      }
      else if (t == NHistYear+1){
        DF <- 1
      }
      else {
        DF <- DF/(1 + popgr + input$eta*incgr + input$PRTP/100)
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
    
    SCCv <- matrix(SCC,NYear,1)
    #                     1     2          3       4         5           6          7       8          9        10        11      12      13      14        15      16                
    resultsout <- cbind(Year, CO2emit, CO2conc, TempAtm, Temperature, population, income, energyint, carbonint, pophist, inchist, eihist, cihist, co2hist, ImpactP, SCCv)
  })
  
    resultsimp <- reactive({
      impdata <- read.csv(file="totalimpact.csv", header=FALSE, sep=",")
      temp <- impdata[1:27,1]
      impobs <- impdata[1:27,2]
      if (input$radio == "1") {
        imp1 <- 0
        imp13 <- -1.0/2.5^.13
        imp2 <- 0
        imp6 <- 0
      }
      else if (input$radio == "2") {
        imp1 <- 0
        imp13 <- 0
        imp2 <- -1.0/3/3
        imp6 <- 0
      }
      else if (input$radio == "3") {
        imp1 <- 0
        imp13 <- 0
        imp2 <- -2.5/3/3
        imp6 <- 0
      }
      else if (input$radio == "4"){
        imp1 <- 2.46
        imp13 <- 0
        imp2 <- -1.11
        imp6 <- 0
      }
      else if (input$radio == "5"){
        imp1 <- -0.12
        imp13 <- 0
        imp2 <- -0.16
        imp6 <- 0
      }
      else if (input$radio == "6"){
        imp1 <- 0
        imp13 <- 0
        imp2 <- -100/20.46/20.46
        imp6 <- -100/6.081^6.754
      }
      else {
        imp1 <- 0
        imp13 <- 0
        imp2 <- 0
        imp6 <- 0
      }
      impmod <- matrix(NA,27,1)
      impmod1 <- matrix(NA,27,1)
      tempgrid <- matrix(NA,27,1)
      tempgrid[1] <- -1
      for (i in 2:27){
        tempgrid[i] <- tempgrid[i-1] + 7/26
      }
      
      for (i in 1:27){
        impmod[i] <- imp1*tempgrid[i] + imp13*abs(tempgrid[i])^1.3 + imp2*tempgrid[i]^2 + imp6*abs(tempgrid[i])^6.754
        impmod1[i] <- imp1*temp[i] + imp13*abs(tempgrid[i])^1.3 + imp2*temp[i]^2 + imp6*abs(temp[i])^6.754
      }
      resultsimp <- cbind(temp,impobs,tempgrid,impmod)  
    })
    
    resultsco2 <- reactive({
      resultsco2 <- read.csv(file="co2concentration.csv", header=FALSE, sep=",")    
    })
    
    output$plotI <- renderPlot({
      plot(resultsimp()[1:27,1],resultsimp()[1:27,2], xlab = 'degree Celsius', ylab = 'welfare equivalent income loss')
      lines(resultsimp()[1:27,3],resultsimp()[1:27,4], type="l", lty="solid")
      title(paste("Impact of climate change"))
    })

    output$plotT <- renderPlot({
      paris <- rep.int(2, NYear)
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,4], type = "l", xlab = 'year', ylab = 'degree Celsius',ylim = c(-0.5,6.0))
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,5]+0.3, type="p")
      lines(resultsout()[1:NYear,1],paris, type="l", lty="dashed")
      title(paste("Temperature"))
    })
    
    output$plotE <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,2]/1000, type = "l", xlab = 'year', ylab = 'billion tonnes of carbon per year')
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,14]/1000, type="p")
      title(paste("Carbon dioxide emissions"))
    })
    
    output$plotP <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,6]/1000000000, type = "l", xlab = 'year', ylab = 'billion people')
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,10]/1000000000, type="p")
      title(paste("Population"))
    })
    
    output$plotY <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,7], type = "l", xlab = 'year', ylab = 'dollar (2010) per person per year')
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,11], type="p")
      title(paste("Per capita income"))
    })
    
    output$plotEI <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,8]*1000, type = "l", xlab = 'year', ylab = 'gram oil equivalent per dollar', ylim = c(0,300))
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,12]*1000, type="p")
      title(paste("Energy intensity"))
    })

    output$plotCI <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,9], type = "l", xlab = 'year', ylab = 'gram carbon per gram oil equivalent')
      lines(resultsout()[1:NHistYear,1],resultsout()[1:NHistYear,13], type="p")
      title(paste("Carbon intensity"))
    })
    
    output$plotC <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,3], type = "l", xlab = 'year', ylab = 'parts per million')
      lines(resultsco2()[,1],resultsco2()[,2], type="p")
      title(paste("Carbon dioxide concentration"))
    })
    
    output$plotD <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,15], type = "l", xlab = 'year', ylab = 'percent income')
      title(paste("Impact of climate change"))
    })
    
    output$SCC <- renderText({paste("Social cost of carbon: $",format(resultsout()[NYear,16]/1000000000,digits=4),"/tC")})
}

# Run the application 
shinyApp(ui = ui, server = server)