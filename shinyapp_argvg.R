setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets directory to current filepath
source("sim_arGamma.R") #load function that simulates an AR gamma process
library(shiny)
library(tidyverse)
library(RColorBrewer)
library(moments)

####ggplot theme####
custom_theme <- function() {
  palette <- brewer.pal("Greys", n=9)
  color.background = palette[1]
  color.grid.major = palette[3]
  color.axis.text = palette[6]
  color.axis.title = palette[7]
  color.title = palette[9]
  
  theme_bw(base_size=9) +
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=7,color=color.axis.title)) +
    theme(plot.title=element_text(color=color.title, size=10, vjust=1.25)) +
    theme(axis.text.x=element_text(size=7,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=7,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=8,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=8,color=color.axis.title, vjust=1.25)) +
    theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}


ui_argvg <- fluidPage(
  
  # App title ----
  titlePanel("Example of the ARGVG model"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for rho ----
      sliderInput(inputId = "rho_r",
                  label = HTML("&rho;:"),
                  min = 0.70,
                  max = 0.99,
                  value = 0.99),
      
      # Input: Slider for nu ----
      sliderInput(inputId = "nu",
                  label = HTML("&nu;:"),
                  min = 0.01,
                  max = 4,
                  value = 2),
      
      # Input: Slider for n ----
      sliderInput(inputId = "n",
                  label = "n:",
                  min = 100,
                  max = 10000,
                  value = 3000)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      fluidRow(
        column(6,plotOutput(outputId="yPlot", width="300px",height="300px")),
        column(6,plotOutput(outputId="ysqPlot", width="300px",height="300px")),
        column(6,plotOutput(outputId="rPlot", width="300px",height="300px")),
        column(6,plotOutput(outputId="yHist", width="300px",height="300px"))
      )
    )
  )
)

server_argvg <- function(input, output) {
  
  
  x    <- reactive({
    rnorm(input$n)
    })
  r <- reactive({
    sim_ARgamma(input$nu, input$rho_r, input$n)
    })
  
  df <- reactive({
    data.frame(r = r(), y = sqrt(r())*x(), ysq = r()*x()^2)
    })
  
  

  output$rPlot <- renderPlot({

    ggplot(df()) +
      geom_line(aes(y = r, x = seq(1:length(df()$r)))) +
      labs(x = "t", y = "R", title = "R") +
      custom_theme()
  })
  
  output$yPlot <- renderPlot({
    ggplot(df()) +
      geom_line(aes(y = y, x = seq(1:length(df()$r)))) +
      labs(x = "t", y = "Y", title = "Y") +
      custom_theme()
  })
  
  output$ysqPlot <- renderPlot({
    ggplot(df()) +
      geom_line(aes(y = ysq, x = seq(1:length(df()$r)))) +
      labs(x = "t", y = "Y squared", title = "Y squared") +
      custom_theme()
  })
  
  output$yHist <- renderPlot({
    ggplot(df(), aes(x=y, y=..density..)) +
      geom_histogram(bins = 50) +
      custom_theme() +
      labs(x = "Y", y = "Frequency", title = "Histogram of Y")
  })
  
}

shinyApp(ui = ui_argvg, server = server_argvg)
