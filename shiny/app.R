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
    column(6,numericInput("VUR", "Cov(YU,YR) (default assumes YR is efficient, and sets to sigmaR^2):", value = 0.09^2, step = 0.1)),
           column(6, actionButton("compute", "Adapt"))
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
                 style = "margin-top: 10px;",
                 condition = "input.compute > 0",
                 h4("Explanation for Adaptation"),
                 tags$ul(
                   tags$li(textOutput("corr_output")),
                   tags$li(textOutput("st_output")),
                   tags$li(textOutput("estimate_output"))
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
  
  regret_st <- round(max(risk_function_st/risk_oracle),2)-1
  
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
  # Show the values in an HTML table ----
  output$values <- renderTable({
    data.frame(
      Estimator = c("Estimate", "Std Error", "Max Regret", "Threshold"),
      Y_U = c(YU, paste("(",sigmaU,")",sep=""), paste(100*round(1/(1-corr^2)-1,2), "%", sep=""), NA),
      Y_R = c(YR, paste("(",sigmaR,")",sep=""), "∞", NA),
      # Adaptive = c(0.36, NA, "44%", NA),
      Soft_threshold = c(adaptive_st, NA, paste(regret_st*100, "%", sep=""), round(st,2)),
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
      labs(x = TeX("b/$\\sigma_O$"), y = TeX("MSE Relative to $\\sigma^2_U$"), title = NULL) +
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
          "If YR is subject to potential bias, then the line labeled ''oracle'' in the figure below illustrates the smallest possible risk that can be achieved when only a bound on the bias magnitude is known.",
          "Adaptation seeks to minimize the maximum regret, which is the worst-case deviation from this oracle risk function.")
  })
  output$st_output <- renderText({
    paste("Given the relative efficiency of YU and the efficient GMM, the soft threshold that achieves optimal adaptation is ", paste(round(st,2),".",sep=""), "The maximum regret in this case is",paste(regret_st*100, "%.", sep=""))
  })
  
  if (-st < tO & tO < st) {
    output$estimate_output <- renderText({
      paste("Since (YR-YU)/std(YR-YU) = ",tO,
            "does not exceed the threshold, the adaptive soft thresholding estimator maintains to be the efficient GMM estimate, yielding",
            paste(adaptive_st,".",sep=""))
      })
  }
  if (tO <=-st ) {
    output$estimate_output <- renderText({
      paste("Since (YR-YU)/std(YR-YU) = ",tO,
            "is negative and below the threshold, the adaptive soft thresholding estimator translates YU toward the GMM efficient estimate, yielding",
            paste(adaptive_st,".",sep=""))
    })
  }
  if (tO >= st ) {
    output$estimate_output <- renderText({
      paste("Since (YR-YU)/std(YR-YU) = ",tO,
            "is positive and above the threshold, the adaptive soft thresholding estimator translates YU toward the efficient GMM estimate, yielding",
            paste(adaptive_st,".",sep=""))
    })
  }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
