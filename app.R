# Shiny Integrated Assessment Model

library(shiny)

# Define UI for user input
ui <- fluidPage(
  titlePanel("Hit the Paris target: Keep temperature below 2 degrees above "),
  fluidRow(
    column(4,
      sliderInput("ClimSens",
        "Warming per CO2 doubling",
        min = 0.1,
        max = 10,
        value = 3.5),
      sliderInput("PeakYear",
        "Year emissions peak",
        min = 2020,
        max = 2090,
        value = 2050),
      sliderInput("DeclineRate",
        "Rate emissions fall",
        min = -0.05,
        max = 0.02,
        value = 0.00)
    ),
    column(4,
      plotOutput("plotH"),
      plotOutput("plotE")
    ),
    column(4,
      plotOutput("plotC"),
      plotOutput("plotT")
    )
  )
)

# Define server to generate output
server <- function(input, output) {
  
  resultsout <- reactive({
    
    source("initialize.R")
    
    RFpar <- input$ClimSens/l2/log(2)
    pkyr <- input$PeakYear-2008 
    decline <- -input$DeclineRate
    histgrowth <- (CO2emit[NHistYear]/CO2emit[NHistYear-10])^0.1 - 1 
    
    for (t in 2:NYear){
      if (t > NHistYear){
        Year[t] <- Year[t-1]+1
        if (t < NHistYear + pkyr){
          growth <- histgrowth * (1 - (NHistYear+pkyr-t)/(pkyr))
          CO2emit[t] <- CO2emit[t-1]*(1+growth)
        }
        else {
          CO2emit[t] <- CO2emit[t-1]*(1-decline)
        }
      }
      for (b in 1:NCO2Box){
        CO2box[t,b] <- CO2box[t-1,b]*(1-BoxLife[b]) + BoxShare[b]*ConvF*CO2emit[t-1]
      }
      CO2conc[t] <- sum(CO2box[t,1:NCO2Box])
      RadForc[t] <- RFpar*log(CO2conc[t]/CO2pre)
      TempAtm[t] <- TempAtm[t-1] + l1*(l2*RadForc[t]-TempAtm[t-1]) + l3*(TempOc[t-1]-TempAtm[t-1])
      TempOc[t] <- TempOc[t-1] + l4*(TempAtm[t-1]-TempOc[t-1])
    }
    resultsout <- cbind(Year, CO2emit, CO2conc, TempAtm, Temperature)
  })

    output$plotH <- renderPlot({
      resid <- resultsout()[1:NYear,4] - resultsout()[1:NYear,5]
      RSS <- 0
      for (t in 101:NHistYear){
        RSS <- RSS + (resid[t]-0.3)^2
      }
      plot(resultsout()[101:NHistYear,1],resultsout()[101:NHistYear,4], type = "l", xlab = 'year', ylab = 'degree Celsius',ylim = c(-0.5,1.5))
      lines(resultsout()[101:NHistYear,1],resultsout()[101:NHistYear,5]+0.3, type="p")
      #legend(1750, 1.5, legend=c("Observed", "Fitted"), lty= c(0,1), pch = c(1,NA))
      title(paste("Observed v modelled temperature, RSS = ", format(RSS,digits=4)))
    })
    
    output$plotE <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,2], type = "l", xlab = 'year', ylab = 'GtC')
      title(paste("Carbon dioxide emissions in 2100:",format(resultsout()[NYear,2],digits=5)))
    })
    output$plotC <- renderPlot({
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,3], type = "l", xlab = 'year', ylab = 'parts per million')
      title(paste("Carbon dioxide concentration in 2100:",format(resultsout()[NYear,3],digits=4)))
    })
    output$plotT <- renderPlot({
      paris <- rep.int(2, NYear)
      plot(resultsout()[1:NYear,1],resultsout()[1:NYear,4], type = "l", xlab = 'year', ylab = 'degree Celsius', ,ylim = c(0.0,8.0))
      lines(resultsout()[1:NYear,1],paris, type="l", lty="dashed")
      title(paste("Temperature in 2100:",format(resultsout()[NYear,4],digits=3)))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)