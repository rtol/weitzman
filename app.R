#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

# Define UI for application that draws a histogram
#ui <- fluidPage(
ui <- navbarPage("Weitzman",
tabPanel("Abatement costs",
        # Application title
    titlePanel("Prices v quantities"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("shift",
                        "Bias in abatement costs",
                        min = -10,
                        max = 10,
                        value = 0,
                        step = 2),
            helpText("The regulator may not know the true cost of emission reduction. Shift the slider to the right if the regulator thinks that abatement is more expensive than it really is, and to left if she thinks it cheaper."),
            sliderInput("steep",
                        "Steepness of benefit curve",
                        min = 0.4,
                        max = 1.6,
                        value = 1,
                        step = 0.2),
            helpText("Both the marginal costs of emission reduction and its marginal benefits vary with emissions. Shift the slider to the right if the marginal benefits are more sensitive to emissions, and to the left if marginal costs are."),
            tags$a(href="https://scholar.harvard.edu/weitzman/files/prices_vs_quantities.pdf", "Read Weitzman's seminal paper."),
            tags$br(),
            tags$a(href="https://www.amazon.co.uk/Climate-Economics-Economic-Analysis-Change/dp/1786435098/", "Any good textbook of environmental economics has an exposition on the Weitzman Theorem."),
            tags$br(),
            tags$a(href="https://github.com/rtol/weitzman", "Code on GitHub")
        ),
        
        mainPanel(
            plotOutput("distPlot"),
            h4(textOutput("reg")),
            h4(textOutput("welf")),
            h4(textOutput("wman"))
        )
    )
),
tabPanel("Damage costs",
    titlePanel("Prices v quantities"),
    
    sidebarLayout(
        sidebarPanel(
            sliderInput("shift2",
                        "Bias in damage costs",
                        min = -10,
                        max = 10,
                        value = 0,
                        step = 2),
            helpText("The regulator may not know the true cost of environmental. Shift the slider to the right if the regulator thinks that abatement is more expensive than it really is, and to left if she thinks it cheaper."),
            sliderInput("steep2",
                        "Steepness of benefit curve",
                        min = 0.4,
                        max = 1.6,
                        value = 1,
                        step = 0.2),
            helpText("Both the marginal costs of emission reduction and its marginal benefits vary with emissions. Shift the slider to the right if the marginal benefits are more sensitive to emissions, and to the left if marginal costs are."),
            tags$a(href="https://scholar.harvard.edu/weitzman/files/prices_vs_quantities.pdf", "Read Weitzman's seminal paper."),
            tags$br(),
            tags$a(href="https://www.amazon.co.uk/Climate-Economics-Economic-Analysis-Change/dp/1786435098/", "Any good textbook of environmental economics has an exposition on the Weitzman Theorem."),
            tags$br(),
            tags$a(href="https://github.com/rtol/weitzman", "Code on GitHub")
        ),
        
        mainPanel(
            plotOutput("distPlot2"),
            h4(textOutput("reg2")),
            h4(textOutput("welf2")),
            h4(textOutput("wman2"))
        )
    )
)
)

# Define server logic
server <- function(input, output) {
    
    results <- reactive({ #mistakes with abatement costs
    
        N <- 101
        step <- 0.1
        emission <- matrix(0,N,1)
        for (i in 2:N){
            emission[i]=emission[i-1]+step
        }
        margcost <- matrix(0,N,1)
        margcost2 <- matrix(0,N,1)
        margben <- matrix(0,N,1)
    
        mc0 <- 20
        mcp <- mc0 + input$shift
        mc1 <- -2
        mb1 <- -mc1*input$steep
        mb0 <- 1
        
        margcost <- mc0 + mc1*emission
        margcost2 <- mcp + mc1*emission
        margben <- mb0 + mb1*emission
        
        optq <- matrix((mc0-mb0)/(mb1-mc1),N,1)                                 # q*
        optp <- matrix(mc0+(mc0-mb0)/(mb1-mc1)*mc1,N,1)                         # p*
        subqq <- matrix((mcp-mb0)/(mb1-mc1),N,1)                                # q'
        subpp <- matrix(mcp+(mcp-mb0)/(mb1-mc1)*mc1,N,1)                        # p'
        subqp <- matrix((mcp+(mcp-mb0)/(mb1-mc1)*mc1-mc0)/mc1,N,1)              # q'''
        subpq <- matrix(mc0+mc1*(mcp-mb0)/(mb1-mc1),N,1)                        # p''
        subpqp <- matrix(mb0+mb1*(mcp+(mcp-mb0)/(mb1-mc1)*mc1-mc0)/mc1,N,1)     # p'''
        
        results <- cbind(emission,margcost,margcost2,margben,optq,optp,subqq,subpp,subqp,subpq,subpqp)
    })

    results2 <- reactive({ #mistakes with damage costs
        
        N <- 101
        step <- 0.1
        emission <- matrix(0,N,1)
        for (i in 2:N){
            emission[i]=emission[i-1]+step
        }
        margcost <- matrix(0,N,1)
        margben <- matrix(0,N,1)
        margben2 <- matrix(0,N,1)
        
        mc0 <- 20
        mc1 <- -2
        mb1 <- -mc1*input$steep2
        mb0 <- 1
        mbp <- mb0 + input$shift2
        
        margcost <- mc0 + mc1*emission
        margben <- mb0 + mb1*emission
        margben2 <- mbp + mb1*emission
        
        optq <- matrix((mc0-mb0)/(mb1-mc1),N,1)                                 # q*
        optp <- matrix(mc0+(mc0-mb0)/(mb1-mc1)*mc1,N,1)                         # p*
        
        subq <- matrix((mc0-mbp)/(mb1-mc1),N,1)                                 # q'
        subp <- matrix(mc0+(mc0-mbp)/(mb1-mc1)*mc1,N,1)                         # p'
        
        subpp <- matrix(mbp+mb1*(mc0-mb0)/(mb1-mc1),N,1)                        # p''
        
        results2 <- cbind(emission,margcost,margben,margben2,optq,optp,subq,subp,subpp)
    })
    
    output$distPlot <- renderPlot({ #mistakes with abatement
        df = data.frame(emission=results()[1:101,1],cost=results()[1:101,2],pcost=results()[1:101,3],benefit=results()[1:101,4])

        harberger <- data.frame(
            loss = c("tax", "tax", "tax", "cap", "cap", "cap"),
            q = c(results()[1,5],results()[1,7],results()[1,7],results()[1,5],results()[1,9],results()[1,9]),  #q* q' q' q* q''' q'''
            p = c(results()[1,6],results()[1,8],results()[1,10],results()[1,6],results()[1,8],results()[1,11]) #p* p' p'' p* p' p''
        )
        
        colover <- c("brown","green")
        colunder <- c("green","brown")
        
        baseplot <- ggplot(data = df,aes(x=emission)) + ylim(-10,30) +
            geom_line(aes(y=cost), color = "brown", size=1) +
            geom_line(aes(y=pcost), color = "red", size=1) +
            geom_line(aes(y=benefit), color = "green", size=1) +
            geom_segment(aes(x=0,y=-10,xend=0,yend=30)) +
            geom_segment(aes(x=0,y=0,xend=10,yend=0)) +
            geom_segment(aes(x=results()[1,5],y=0,xend=results()[1,5],yend=results()[1,6])) +
            geom_segment(aes(x=0,y=results()[1,6],xend=results()[1,5],yend=results()[1,6])) +
            annotate(geom="text", x=results()[1,5], y=-1, label="q*", color="black") +
            annotate(geom="text", x=-0.2, y=results()[1,6], label="p*", color="black")
            
        if (input$shift > 0){
            weitzplot <- baseplot + 
                geom_segment(aes(x=results()[1,7],y=0,xend=results()[1,7],yend=results()[1,8]),linetype = "dashed") +
                geom_segment(aes(x=0,y=results()[1,8],xend=results()[1,7],yend=results()[1,8]),linetype = "dashed") +
                annotate(geom="text", x=results()[1,7], y=-1, label="q'", color="black") +
                annotate(geom="text", x=-0.2, y=results()[1,8], label="p'", color="black") +
                geom_segment(aes(x=results()[1,9],y=0,xend=results()[1,9],yend=results()[1,8]),linetype = "dotted") +
                geom_segment(aes(x=0,y=results()[1,10],xend=results()[1,7],yend=results()[1,10]),linetype = "dotted") +
                annotate(geom="text", x=results()[1,9], y=-1, label="q'''", color="black") +
                geom_polygon(data = harberger, mapping = aes(x = q, y = p, group = loss, fill = loss)) +
                scale_fill_manual(values = colover)
            if (input$steep == 1){
                weitzplot +
                    annotate(geom="text", x=-0.2, y=results()[1,10], label="p'''", color="black")
            }
            else {
                weitzplot +
                    geom_segment(aes(x=0,y=results()[1,11],xend=results()[1,9],yend=results()[1,11]),linetype = "dotted")  +
                    annotate(geom="text", x=-0.2, y=results()[1,10], label="p''", color="black") +
                    annotate(geom="text", x=-0.2, y=results()[1,11], label="p'''", color="black")
            }
        }
        else if (input$shift < 0){
            color <- c("green")
            weitzplot <- baseplot + 
                geom_segment(aes(x=results()[1,7],y=0,xend=results()[1,7],yend=results()[1,8]),linetype = "dashed") +
                geom_segment(aes(x=0,y=results()[1,8],xend=results()[1,7],yend=results()[1,8]),linetype = "dashed") +
                annotate(geom="text", x=results()[1,7], y=-1, label="q'", color="black") +
                annotate(geom="text", x=-0.2, y=results()[1,8], label="p'", color="black") +
                geom_segment(aes(x=results()[1,9],y=0,xend=results()[1,9],yend=results()[1,8]),linetype = "dotted") +
                geom_segment(aes(x=0,y=results()[1,11],xend=results()[1,9],yend=results()[1,11]),linetype = "dotted") +
                annotate(geom="text", x=results()[1,9], y=-1, label="q'''", color="black") +
                geom_polygon(data = harberger, mapping = aes(x = q, y = p, group = loss, fill = loss)) +
                scale_fill_manual(values = colunder)
            if (input$steep == 1){
                weitzplot  +
                    annotate(geom="text", x=-0.2, y=results()[1,10], label="p'''", color="black")
            }
            else {
                weitzplot +
                    geom_segment(aes(x=0,y=results()[1,10],xend=results()[1,7],yend=results()[1,10]),linetype = "dotted")  +
                    annotate(geom="text", x=-0.2, y=results()[1,10], label="p''", color="black") +
                    annotate(geom="text", x=-0.2, y=results()[1,11], label="p'''", color="black")
            }
        }
        else{
            baseplot
        }
        
    })

    output$distPlot2 <- renderPlot({ #mistakes with damage
        df = data.frame(emission=results2()[1:101,1],cost=results2()[1:101,2],benefit=results2()[1:101,3],pben=results2()[1:101,4])
        
        harberger <- data.frame(
            loss = c("tax or cap", "tax or cap", "tax or cap"),
            q = c(results2()[1,5],results2()[1,7],results2()[1,5]),  #q* q' q*
            p = c(results2()[1,6],results2()[1,8],results2()[1,9])   #p* p' p''
        )
        
        col <- c("red")
        
        baseplot <- ggplot(data = df,aes(x=emission)) + ylim(-10,30) +
            geom_line(aes(y=cost), color = "brown", size=1) +
            geom_line(aes(y=benefit), color = "green", size=1) +
            geom_line(aes(y=pben), color = "red", size=1) +
            geom_segment(aes(x=0,y=-10,xend=0,yend=30)) +
            geom_segment(aes(x=0,y=0,xend=10,yend=0)) +
            geom_segment(aes(x=results2()[1,5],y=0,xend=results2()[1,5],yend=results2()[1,6])) +
            geom_segment(aes(x=0,y=results2()[1,6],xend=results2()[1,5],yend=results2()[1,6])) +
            annotate(geom="text", x=results2()[1,5], y=-1, label="q*", color="black") +
            annotate(geom="text", x=-0.2, y=results2()[1,6], label="p*", color="black")
        
        if (input$shift2 != 0){
            weitzplot <- baseplot + 
                geom_segment(aes(x=results2()[1,7],y=0,xend=results2()[1,7],yend=results2()[1,8]),linetype = "dashed") +
                geom_segment(aes(x=0,y=results2()[1,8],xend=results2()[1,7],yend=results2()[1,8]),linetype = "dashed") +
                annotate(geom="text", x=results2()[1,7], y=-1, label="q'", color="black") +
                annotate(geom="text", x=-0.2, y=results2()[1,8], label="p'", color="black") +
                geom_segment(aes(x=0,y=results2()[1,9],xend=results2()[1,5],yend=results2()[1,9]),linetype = "dotted") +
                annotate(geom="text", x=-0.2, y=results2()[1,9], label="p''", color="black") +
                geom_polygon(data = harberger, mapping = aes(x = q, y = p, group = loss, fill = loss)) +
                scale_fill_manual(values = col)
            weitzplot
            
        }
        
        else{
            baseplot
        }
        
    })
        
    output$reg <- renderText({ #mistake with abatement
        if (input$shift == 0) {
            "If the marginal costs of emission reduction and its marginal benefits are known, it does not matter whether the regulator picks a price instrument (taxes) or quantity instrument (cap-and-trade). A tax at level p* leads to emissions q*, and a cap at level q* implies a price p*."
        }
        else if (input$shift > 0){
            "However, if the regulator thinks that emissions reduction is more expensive than it really is, price and quantity instruments are no longer equivalent. In this case, the regulator would set the cap at q', that is, allow too high emissions, and the tax at p', that is, impose too high a tax."
        }
        else {
            "However, if the regulator thinks that emissions reduction is less expensive than it really is, price and quantity instruments are no longer equivalent. In this case, the regulator would set the cap at q', that is, reduce emissions too far, and the tax at p', that is, impose too low a tax."
        }
    })
    
    output$reg2 <- renderText({ #mistake with damage
        if (input$shift2 == 0) {
            "If the marginal costs of emission reduction and its marginal benefits are known, it does not matter whether the regulator picks a price instrument (taxes) or quantity instrument (cap-and-trade). A tax at level p* leads to emissions q*, and a cap at level q* implies a price p*."
        }
        else if (input$shift2 > 0){
            "If the regulator thinks that environmental damage is greater than it really is, she would set the cap at q', that is, reduce emissions too far, and the tax at p', that is, impose too high a tax. Both price and quantity instruments lead to overregulation."
        }
        else {
            "If the regulator thinks that environmental damage is smaller than it really is, she regulator would set the cap at q', that is, allow too high emissions, and the tax at p', that is, impose a tax that is too low. Both price and quantity instruments lead to underregulation."
        }
    })
    
    output$welf <- renderText({  #mistake with abatement
        if (input$shift == 0) {
            "There is no welfare loss, as regulation is set at the level where marginal costs equal marginal benefits."
        }
        else {
            "Overregulation leads to a welfare gain for those who care about the environment, but a welfare loss for the polluters. Underregulation implies a welfare gain for the polluters, and a welfare loss for the environment. The total welfare loss is shown in the diagram, in green if the environment gains more than the polluters lose, in brown otherwise. The legend indicates whether the welfare loss is due to an inappropriate tax or cap."
        }
    })

    output$welf2 <- renderText({  #mistake with damage
        if (input$shift2 == 0) {
            "There is no welfare loss, as regulation is set at the level where marginal costs equal marginal benefits."
        }
        else if (input$shift2 > 0){
            "Overregulation leads to a welfare gain for those who care about the environment, but a welfare loss for the polluters. The total welfare loss is shown in the diagram."
        }
        else {
            "Underregulation implies a welfare gain for the polluters, and a welfare loss for the environment. The total welfare loss is shown in the diagram."
        }
    })
        
    output$wman <- renderText({ #mistake with abatement
        if (input$shift == 0) {
            "The relative steepness of the marginal cost and benefit curves is irrelevant."
        }
        else {
            if (input$steep == 1){
                "If the marginal benefit curve is as steep as the marginal cost curve, the welfare loss due to the wrong cap equals the welfare loss due to the wrong tax, although the distribution of welfare is different."
            }
            else if (input$steep >1){
                "If the marginal benefit curve steeper than the marginal cost curve, the welfare loss due to the wrong cap is smaller than the welfare loss due to the wrong tax. A wise regulator would impose a cap."
            }
            else {
                "If the marginal benefit curve shallower than the marginal cost curve, the welfare loss due to the wrong cap is larger than the welfare loss due to the wrong tax. A wise regulator would impose a tax."
            }
        }
    })
    
    output$wman2 <- renderText({ #mistake with damage
        "The relative steepness of the marginal cost and benefit curves is irrelevant."
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
