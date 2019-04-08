# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
#library(ggplot2)
#library(Matrix)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Choose the climate sensitivity to minimize RSS"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("ClimSens",
                  "Warming per CO2 doubling",
                  min = 0.1,
                  max = 10,
                  value = 2.5),
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
    
    # Show a plot of the generated distribution
    mainPanel("Model output",
              plotOutput("tempPlot")
              #plotOutput("output$plotgraph1")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$tempPlot <- renderPlot({
    
    source("initialize.R")
    
    RFpar <- input$ClimSens/l2/log(2)
    pkyr <- input$PeakYear-2008 
    decline <- -input$DeclineRate
    
    for (t in 2:NYear){
      if (t > NHistYear){
        Year[t] <- Year[t-1]+1
        if (t < NHistYear + pkyr){
          growth <- 0.02 * (1 - (NHistYear+pkyr-t)/(pkyr))
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
    resid <- TempAtm - Temperature
    RSS <- 0
    for (t in 101:NHistYear){
      RSS <- RSS + (TempAtm[t]-Temperature[t]-0.3)^2
    }
    attach(mtcars)
    par(mfrow=c(4,1))
    plot(Year[101:NHistYear],TempAtm[101:NHistYear], type = "l", xlab = 'year', ylab = 'degree Celsius', ylim = c(-0.6,1.5))
    lines(Year[101:NHistYear],Temperature[101:NHistYear]+0.3, type="p")
    legend(1750, 1.5, legend=c("Observed", "Fitted"), lty= c(0,1), pch = c(1,NA))
    title(paste("Observed v modelled temperature, RSS = ", format(RSS,digits=3)))
    plot(Year[1:NYear],CO2emit[1:NYear], type = "l", xlab = 'year', ylab = 'GtC')
    title(paste("Carbon dioxide emissions in 2100:",format(CO2emit[NYear],digits=5)))
    plot(Year,CO2conc, type = "l", xlab = 'year', ylab = 'parts per million')
    title(paste("Carbon dioxide concentration in 2100:",format(CO2conc[NYear],digits=4)))
    plot(Year,TempAtm, type = "l", xlab = 'year', ylab = 'degree Celsius')
    title(paste("Temperature in 2100:",format(TempAtm[NYear],digits=3)))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)