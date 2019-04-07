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
                     "Value:",
                     min = 0.1,
                     max = 10,
                     value = 3.5)#,
         #sliderInput("OceanSpeed",
          #           "Value:",
          #           min = 100,
           #          max = 300,
            #         value = 176)
      ),
      
      # Show a plot of the generated distribution
      mainPanel("Observed and modelled temperature",
         plotOutput("tempPlot")
         #plotOutput("output$plotgraph1")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
     output$tempPlot <- renderPlot({
      #source("run.SIAM.R(input$ClimSens)")
      
      source("initialize.R")
      
  #pt1 <- reactive({
      RFpar <- input$ClimSens/l2/log(2)
      #l4 <- 1/input$OceanSpeed
      
      for (t in 2:NHistYear){
        for (b in 1:NCO2Box){
          CO2box[t,b] = CO2box[t-1,b]*(1-BoxLife[b]) + BoxShare[b]*ConvF*CO2emit[t-1,2]
        }
        CO2conc[t] <- sum(CO2box[t,1:NCO2Box])
        RadForc[t] <- RFpar*log(CO2conc[t]/CO2pre)
        TempAtm[t] <- TempAtm[t-1] + l1*(l2*RadForc[t]-TempAtm[t-1]) + l3*(TempOc[t-1]-TempAtm[t-1])
        TempOc[t] <- TempOc[t-1] + l4*(TempAtm[t-1]-TempOc[t-1])
      }
      resid <- TempAtm - Temperature
      RSS <- 0
      for (t in 101:NHistYear){
        RSS <- RSS + (TempAtm[t]-Temperature[t,2]-0.3)^2
      }
      plot(CO2emit[1:NHistYear,1],TempAtm, type = "l", xlab = 'year', ylab = 'annual mean surface air temperature', ylim = c(-0.6,1.5))
      lines(Temperature[1:NHistYear,1],Temperature[1:NHistYear,2]+0.3, type="p")
      legend(1750, 1.5, legend=c("Observed", "Fitted"), lty= c(0,1), pch = c(1,NA))
      title(paste("RSS = ", RSS))
      
  #})
      #pt2 <- plot(CO2emit[1:NHistYear,1],CO2emit[1:NHistYear,2])
      
      #output$plotgraph1 = renderPlot({pt1})
      #output$plotgraph2 = renderPlot({pt2})
     
      })
}

# Run the application 
shinyApp(ui = ui, server = server)

