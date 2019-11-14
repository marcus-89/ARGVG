setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets directory to current filepath
source("sim_arGamma.R") #load function that simulates an AR gamma process
library(shiny)
library(tidyverse)
library(RColorBrewer)
library(moments)
library(cowplot)


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


ui_est <- fluidPage(
  
  # App title ----
  titlePanel("Parameter Estimation"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for rho ----
      sliderInput(inputId = "rho_r",
                  label = HTML("true &rho;:"),
                  min = 0.70,
                  max = 0.99,
                  value = 0.99),
      
      # Input: Slider for nu ----
      sliderInput(inputId = "nu",
                  label = HTML("true &nu;:"),
                  min = 0.01,
                  max = 4,
                  value = 2),
      
      # Input: Slider for n ----
      sliderInput(inputId = "n",
                  label = "n:",
                  min = 100,
                  max = 10000,
                  value = 3000),
      fluidRow(
        column(6,uiOutput("nu_hat_symbol")),
        column(6,uiOutput("rho_hat_symbol"))
      ),
      fluidRow(
        column(2),
        column(4,textOutput("nu_hat")),
        column(2),
        column(4,textOutput("rho_hat"))
      )
    
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      fluidRow(
        column(12,plotOutput(outputId="yDensity", width="600px",height="300px"))
      ),
      hr(),
      fluidRow(
        column(12,plotOutput(outputId = "regression", width="600px", height = "300px"))
      )
    )
  )
)


server_est <- function(input, output) {
  
  
  x    <- reactive({
    rnorm(input$n)
  })
  r <- reactive({
    sim_ARgamma(input$nu, input$rho_r, input$n)
  })
  
  df <- reactive({
    data.frame(r = r(), y = sqrt(r())*x(), ysq = r()*x()^2)
  })
  
  nu_hat = reactive({
    (moments::kurtosis(df()$y)-3)/3
  })

  dens_norm_df = reactive({
    data.frame(y = dnorm(seq(min(df()$y), max(df()$y), 0.01), mean = mean(df()$y), sd = sd(df()$y)),
               x = seq(min(df()$y), max(df()$y), 0.01))
  })
  
  rho_ysq_hat_full = reactive({
    acf(df()$ysq, lag.max = length(df()$ysq), plot=F)$acf[-1]
  })
  
  req_lags = reactive({
    which(rho_ysq_hat_full() < 2/sqrt(length(df()$y)))[1]
  })
  
  rho_ysq_hat = reactive({
    rho_ysq_hat_full()[1:(req_lags()-1)]
  })
  
  rho_rk_hat = reactive({
    rho_ysq_hat()*(2+3*nu_hat())/nu_hat()
  })
  
  linmod = reactive({
    lm(log(rho_rk_hat())~seq(1:length(rho_rk_hat())))
  })
  
  rho_r_hat = reactive({
    exp(linmod()$coefficients[2])
  })
  
  log_rho_r_hat = reactive({
    linmod()$coefficients[2]
  })
  
  rho_est_df = reactive({
    data.frame(hat = rho_rk_hat(), loghat = log(rho_rk_hat()), lag = 1:length(rho_rk_hat()))
  })
  
  fitted_df = reactive({
    data.frame(fit = exp(linmod()$fitted.values), lag = 1:length(rho_rk_hat()))
  })
  
  gg1 <- reactive({
    ggplot(rho_est_df()) +
      geom_point(aes(x=lag, y=hat)) +
      geom_line(data = fitted_df(), aes(y=fit, x=lag), col="red", size=1.1) +
      custom_theme() +
      labs(x="Lag", title = expression(paste("The exponential decay in ", hat(rho)[r]^"k")), 
           y="Estimated Autocorrelation") +
      annotate("text", label=paste0("hat(rho)[r]==~exp(hat(beta))==", round(rho_r_hat(), 3)),
               x = 0.3*length(rho_rk_hat()), y = 0.2*max(rho_rk_hat()), col="red", size = 4.5, parse=T)
  })
  
  gg2 <- reactive({
    ggplot(rho_est_df()) +
      geom_point(aes(x=lag, y=loghat)) +
      geom_abline(intercept=linmod()$coefficients[1], slope=linmod()$coefficients[2],col="red", size=1.1) +
      custom_theme() +
      labs(x="Lag", y = expression(paste("log(",hat(rho)[r]^k, ")")), 
           title = expression(paste("The linear relationship between lag and ", log(hat(rho)[r]^"k")))) +
      annotate("text", label=paste0("hat(beta)==", round(log(rho_r_hat()), 3)),
               x = 0.3*length(rho_rk_hat()), y = 0.8*min(log(rho_rk_hat())), col="red", size = 4.5, parse=T)
  })
  
  output$regression <- renderPlot({
    plot_grid(gg1(),gg2())
  })
  
  output$yDensity <- renderPlot({
    ggplot(df(), aes(x=y, y=..density..)) +
      geom_histogram(bins = 50, fill = "coral3", alpha = 0.4) +
      geom_density(col="red", size=1.05) +
      geom_line(data = dens_norm_df(), aes(x=x, y=y), size=1.05) +
      custom_theme() +
      labs(x = "Y", y = "Frequency", title = 
             expression(paste("Density of ", Y[t], ". Its kurtosis is used to estimate ",
                              hat(nu)))) +
      annotate("text", label=paste("Kurtosis = ", 
                                   round(moments::kurtosis(df()$y),2)),
               x = 0.8*max(df()$y), y = max(dens_norm_df()$y)+0.1, col="red", size = 4.5) +
      annotate("text", label="Kurtosis = 3",
               x = 0.8*max(df()$y), y = max(dens_norm_df()$y), size = 4.5)
  })
  
  output$nu_hat <- renderText({
    paste(round(nu_hat(),2))
  })
  
  output$rho_hat <- renderText({
    paste(round(rho_r_hat(),3))
  })
  
  output$nu_hat_symbol <- renderUI({
    withMathJax(
      helpText('$$\\LARGE\\hat{\\nu}$$')
    )
  })
  
  output$rho_hat_symbol <- renderUI({
    withMathJax(
      helpText('$$\\LARGE\\hat{\\rho_r}$$')
    )
  })
  
}

shinyApp(ui = ui_est, server = server_est)

