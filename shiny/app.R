# Load necessary libraries
library(shiny)
library(bslib)
source("grid.R")
library(ggplot2)
library(tidyverse)
# Define color and linetype preferences
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
my_linetypes <- c('solid', 'dotted', 'dashed', 'dotdash','longdash')
my_linetypes2 <- c('solid', 'solid', 'dashed', 'dotdash','longdash')

library(latex2exp)
# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("Adaptive Soft Threshold for Combining Unrestricted and Restricted Estimates"),
  # Description line below the title
  h4("Shiny App based on the MissAdapt repository."),
  h4("Adapts to bias due to potential misspecification of the restriction via soft-thresholding."),
  h4("For more functionality refer to the vignette on", a("GitHub.", 
  href = "https://github.com/lsun20/MissAdapt",target="_blank", rel="noreferrer noopener")),
  
 
      fluidRow(
        column(6, numericInput("YU", "Unrestricted estimate (YU):", value = 0.43)),
        column(6, numericInput("sigmaU", "Standard error (sigmaU):", value = 0.14))
      ),
      fluidRow(
        column(6, numericInput("YR", "Restricted estimate (YR):", value = 0.26)),
        column(6, numericInput("sigmaR", "Standard error (sigmaR):", value = 0.09))
      ),
  fluidRow(
    column(6,numericInput("VUR", "Cov(YU,YR): If YU and YR are estimated on independent datasets, please enter zero. The default assumes that YU and YR are estimated from the same dataset, with YR being the efficient estimator. In this case, the covariance is equal to sigmaR^2 (Hausman,1978)."
                          , value = 0.09^2, step = 0.1)),
    column(6,uiOutput("sigmaO"))
    ),

  fluidRow(
    column(6, actionButton("compute", "Adapt")),
    column(6, sliderInput("slider",
                          "Please specify an upper bound on the bias of YR in multiples of std(YR-YU). This bound will be used to construct a fixed-length confidence interval (FLCI). Larger bounds yield longer FLCIs with better worst case coverage.",
                          min = 0,
                          max = 9,
                          value = 1))
  ),
    fluidRow(
             div(
               style = "display: flex; justify-content: center; flex-direction: column; align-items: center;",
               conditionalPanel(
                 condition = "input.compute > 0",
                 h4("Summary of the Adaptation Results"),
                 tableOutput("values")
               )
             )
    ),
    fluidRow(
      column(12,
             div(
               conditionalPanel(
                 style = "list-style-type: disc; margin-top: 10px;",
                 condition = "input.compute > 0",
                 h4("Explanation for Adaptation"),
                 tags$ul(
                   tags$li(style = "margin-bottom: 10px;",textOutput("corr_output")),
                   tags$li(style = "margin-bottom: 10px;",textOutput("estimate_output")),
                   tags$li(style = "margin-bottom: 10px;",textOutput("st_output")),
                   tags$li(style = "margin-bottom: 10px;",textOutput("cv_output"))
                 )
               )
             )
        )
      ),
      fluidRow(
        div(
          style = "display: flex; justify-content: center; flex-direction: column; align-items: center;",
          plotOutput('plot', width = "400px", height = "300px")
        )
      )
    
  
)

# Define server logic
server <- function(input, output,session) {
  observe({
    # Update VUR to default value based on sigmaR
    updateNumericInput(session, "VUR", value = input$sigmaR^2)
    
    
  })
  
  output$sigmaO <- renderUI({
    # Initialize an empty HTML string
    output_html <- ""
    
    if(input$sigmaR >= input$sigmaU) {
      output_html <- "<b>YU is no less precise than YR and already optimal. There is no need to adapt. Please do not proceed.</b>"
    } else {
      bias <- round(input$YR - input$YU, 2)
      std_error <- sqrt((input$sigmaR)^2 - 2*input$VUR + (input$sigmaU)^2)
      
      # Construct the HTML for bias and standard error
      output_html <- paste("<b>The bias in YR has a simple estimate YR-YU = </b>", bias, 
                           ".<br><b>The standard error of this estimate is std(YR-YU) = </b>", 
                           round(std_error,2), ".", sep="")
      
      Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
      corr <- (input$VUR - (input$sigmaU)^2) / std_error / input$sigmaU
      
      output_html <- paste(output_html, "<br><b>The correlation coefficient between YU and (YR-YU) is therefore rho =</b>", paste(round(corr,2),",",sep=""),
          "<b>which implies the relative efficiency of YU  is 1-rho^2 =</b>", paste(round(1-corr^2,2),".",sep=""),
          "<br><b>This relative efficiency determines the amount of adaptation regret.</b>")
      # Check conditions and append to the HTML string
      if(abs(corr) > max(Sigma_UO_grid)) {
        output_html <- paste(output_html,  "<br><b>Error:</b> YR is very precise relative to YU and adaptation incurs large regret. Proceed with caution.")
      }
      
      if(abs(corr) < min(Sigma_UO_grid)) {
        output_html <- paste(output_html,  "<br><b>Error:</b> YU is not much less precise than YR. Adaptation won't do much.:)")
      }
    }
    
    # Return the final HTML
    HTML(output_html)
  })
  observeEvent(input$compute, {
    
    YU <- as.numeric(input$YU)
    sigmaU <- as.numeric(input$sigmaU)
    YR <- as.numeric(input$YR)
    sigmaR <- as.numeric(input$sigmaR)
    VUR <- as.numeric(input$VUR)

    VU <- sigmaU^2
    VR <- sigmaR^2
    VO <- VR - 2 * VUR + VU
    VUO <- (VUR - VU)
    corr <- VUO / sqrt(VO) / sqrt(VU)
  
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
   
  st.function <- splinefun(Sigma_UO_grid, st_thresholds_grid, method = "fmm", ties = mean)
  st <- st.function(abs(corr))
  
  YO <- YR - YU
  tO <- YO / sqrt(VO)
  GMM <- YU - VUO / VO * YO; V_GMM <- VU - VUO/VO*VUO;
  adaptive_st <- VUO / sqrt(VO) * ((tO > st) * (tO - st) + (tO < -st) * (tO + st)) + GMM
  decimal <- match(TRUE, round(YU, 1:20) == YU)
  adaptive_st <- round(adaptive_st,decimal)
  tO <- round(tO,decimal)
  GMM <- round(GMM,decimal)
  # Add soft threshold risk functions
  Eb <- function(l) 1 + l^2 + (b_grid^2 - 1 - l^2) * (pnorm(l - b_grid) - pnorm(-l - b_grid)) + (-b_grid - l) * dnorm(l - b_grid) - (l - b_grid) * dnorm(-l - b_grid)
  
  risk_function_st <- corr^2 * Eb(st) + 1 - corr^2 
  risk_oracle <- corr^2 * rho_b_over_sigma + 1 - corr^2
  
  risk_st <- round(max(risk_function_st),2)-1
  regret_st <- round(max(risk_function_st/risk_oracle),2)-1
  regret_nonlinear_function <- splinefun(Sigma_UO_grid, adaptive_nonlinear_regret_grid, method = "fmm", ties = mean)
  regret_nonlinear <- round(regret_nonlinear_function(abs(corr)),2)-1
  # Use simulation to calculate the risk function for the pre-test estimator that switches btw Y_U and Y_R
  # sims <- 100000
  # x <- rnorm(sims, 0, 1)
  # x_b <- outer(x, as.vector(b_grid), "+")
  # Ebsims_ht <- function(l) {
  #   colMeans(((x_b > l) * x_b + ((x_b < l & x_b > -l) * x_b) * (1 + VO/VUO) + (x_b < -l) * x_b
  #             - matrix(rep(b_grid, sims), ncol = length(b_grid), byrow = TRUE))^2)
  # }
  # risk_function_ht_ttest <- corr^2 * Ebsims_ht(1.96) + 1 - corr^2
  
  # Calculate the risk function for the pre-test estimator that switches btw Y_U and GMM
  # Eb_ht <-function(l){
  #   1+(b_grid^2-1)*(pnorm(l-b_grid)-pnorm(-l-b_grid))+(l-b_grid)*dnorm(l-b_grid) - (-l-b_grid)*dnorm(-l-b_grid);
  # } 
  # risk_function_ht_ttest <- corr^2 * Eb_ht(1.96) + 1 - corr^2
  
  # regret_ht <- round(max(risk_function_ht_ttest/risk_oracle),2)-1
  
  ### Calculate the coverage distortion
  Sigma_UO_grid <- tanh(seq(-3, -0.05, 0.05))
  B_FLCI <- input$slider
  
  flci_st_c_function <-  splinefun(Sigma_UO_grid, flci_adaptive_st_cv_grid[B_FLCI+1,], method = "fmm", ties = mean) 
  c <- flci_st_c_function(corr)
  
  Kb <- length(b_grid)
  corr2 <- corr^2
  rej_prob <- c(Kb)
  for (j in 1:Kb) {
    b <- b_grid[j]
    pos = function (t) ( 1-pnorm((c-corr*(t-st-b))/sqrt(1-corr2)) + pnorm((-c-corr*(t-st-b))/sqrt(1-corr2)) ) * dnorm(t-b)
    neg = function (t) ( 1-pnorm((c-corr*(t+st-b))/sqrt(1-corr2)) + pnorm((-c-corr*(t+st-b))/sqrt(1-corr2)) ) * dnorm(t-b)
    zero = function (t) ( 1-pnorm((c-corr*(-b))/sqrt(1-corr2)) + pnorm((-c-corr*(-b))/sqrt(1-corr2)) ) * dnorm(t-b)
    
    rej_prob[j] <- integrate(pos, lower = st, upper = Inf)$value +
      integrate(zero, lower = -st, upper = st)$value + 
      integrate(neg, lower = -Inf, upper = -st)$value
  }
  flci_title <- TeX("1*$\\sigma_O$-FLCI")
  # Show the values in an HTML table ----
  output$values <- renderTable({
    data.frame(
      Estimator = c("Estimate", "Std Error", "Threshold", "Max Regret", "Max Risk", "FLCI" , "Best coverage", "Worst coverage"),
      Y_U = c(YU, paste("(",sigmaU,")",sep=""), NA, paste(100*round(1/(1-corr^2)-1,2), "%", sep=""), "0%",
              paste("[",round(YU-1.96*sigmaU,decimal),",",round(YU+1.96*sigmaU,decimal),"]", sep=""),"95%","95%"),
      Y_R = c(YR, paste("(",sigmaR,")",sep=""), NA, "∞", "∞", NA, NA, NA),
      # Adaptive = c(0.36, NA, "44%", NA),
      Soft_threshold = c(adaptive_st, NA,  round(st,2),paste(regret_st*100, "%", sep=""),paste(risk_st*100, "%", sep=""),
                         paste("[",round(adaptive_st-c*sigmaU,decimal),",",round(adaptive_st+c*sigmaU,decimal),"]", sep=""),
                         paste(100*round(1-min(rej_prob),2), "%", sep=""), paste(100*round(1-max(rej_prob),2), "%", sep="")),
      # Pre_test = c(0.26, NA, paste(regret_ht*100, "%", sep=""), "1.96"),
      stringsAsFactors = FALSE)
  })
  
  # Create a data frame of only adaptive risk
  # Plot the figure
  b_grid_scale <- sign(b_grid)*sqrt(abs(b_grid));
  df <- data.frame(
    b_grid_scale = b_grid_scale,
    a = rep(1, Kb),
    b = risk_function_st,
    # c = risk_function_ht_ttest,
    e = risk_oracle
  )
  
  # Melt the data frame to long format
  df_long <- tidyr::gather(df, key = "series", value = "y", -b_grid_scale)
  
  output$plot <- renderPlot({
    par(mar = c(5, 5, 4, 5), pty = 'm',cex=1.2)  # Adjust margins as needed
    
    # Plotting using ggplot2
    ggplot(df_long, aes(x = b_grid_scale, y = y, linetype = series, color = series)) +
      geom_line(size=1.2) +
      labs(x = TeX("b/std(YR-YU)"), y = TeX("Risk (MSE) Relative to $\\sigma^2_U$"), title = NULL) +
      scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3), labels = c('-9', '-4', '-1', '0', '1', '4', '9')) +
      ylim(c(min(risk_oracle-0.1), max(risk_function_st)+0.1)) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.text = element_text(size = 12),  # Adjust legend text size
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)
      ) +
      scale_color_manual(values = cbp1, name = "",
                         labels = c(TeX('$Y_U$'),TeX('Adaptive soft-threshold'),'Oracle')) +
      scale_linetype_manual(values = my_linetypes2, name = "",
                            labels = c(TeX('$Y_U$'),TeX('Adaptive soft-threshold'),'Oracle')) +
      guides(color = guide_legend(ncol = 3))
  })
  
  output$corr_output <- renderText({
    paste("The correlation coefficient between YU and (YR-YU) is ", paste(round(corr,2),",",sep=""),
          "which implies the efficient GMM estimate when YR is correctly specified is ",GMM,
          " and the efficiency of YU relative to the efficient GMM estimate is ", paste(round(1-corr^2,2),".",sep=""),
          "If YR is subject to potential bias, then the line labeled ''Oracle'' in the figure below illustrates the smallest possible risk that can be achieved when only a bound on the bias magnitude is known.",
          "Adaptation seeks to minimize the maximum regret, which is the worst-case deviation from this oracle risk function.")
  })
  output$st_output <- renderText({
    paste("The risk of this adaptive soft-threshold estimator is shown in the figure to illustrate the large efficiency gain when bias is small.",
          "The maximum regret in this case is",paste(regret_st*100, "%.", sep=""),"The closest performance one can get relative to the oracle is",
          paste(regret_nonlinear*100, "%,", sep=""), "achieved using a fully nonlinear estimator.","Soft-thresholding provides a simple approximation with minimal efficiency loss." )
  })
  output$cv_output <- renderText({
    paste("If the bias in YR is no larger than",B_FLCI,"*std(YR-YU), a critical value of ",  round(c,2) , 
          "(instead of 1.96) centered at the adaptive soft threshold estimate with std = sigmaU ensures 95% coverage. This CI is different from the one centered at YU due to the bias-efficiency trade-off. ",
          "If the bias is unrestricted, this CI can be used with consideration of the best and worst coverage scenarios.",
          "In contrast, the CI centered at YU has 95% coverage, regardless of the bias.")
  })
  if (-st < tO & tO < st) {
    output$estimate_output <- renderText({
      paste("Given the relative efficiency of YU and GMM, the soft threshold that achieves optimal adaptation is ", paste(round(st,2),".",sep=""), 
            "Since (YR-YU)/std(YR-YU) = ",tO,
            "does not exceed the threshold, the adaptive soft thresholding estimator maintains to be the GMM estimate, yielding",
            paste(adaptive_st,".",sep=""))
      })
  }
  if (tO <=-st ) {
    output$estimate_output <- renderText({
      paste("Given the relative efficiency of YU and GMM, the soft threshold that achieves optimal adaptation is ", paste(round(st,2),".",sep=""), 
            "Since (YR-YU)/std(YR-YU) = ",tO,
            "is negative and below the threshold, the adaptive soft thresholding estimator translates YU towards the GMM estimate, yielding",
            paste(adaptive_st,".",sep=""))
    })
  }
  if (tO >= st ) {
    output$estimate_output <- renderText({
      paste("Given the relative efficiency of YU and GMM, the soft threshold that achieves optimal adaptation is ", paste(round(st,2),".",sep=""), 
            "Since (YR-YU)/std(YR-YU) = ",tO,
            "is positive and above the threshold, the adaptive soft thresholding estimator translates YU towards the GMM estimate, yielding",
            paste(adaptive_st,".",sep=""))
    })
  }
  
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
