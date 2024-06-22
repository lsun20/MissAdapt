# Load necessary libraries
library(shiny)

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("Adaptive Soft Threshold for Combining Unrestricted and Restricted estimates"),
  # Description line below the title
  h4("Shiny app based on the MissAdapt function. For more functionality visit: ", 
     a("https://github.com/lsun20/MissAdapt", 
       href = "https://github.com/lsun20/MissAdapt",target="_blank")),
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(6, numericInput("YU", "Unrestricted estimate (YU):", value = 0.0043)),
        column(6, numericInput("sigmaU", "Standard error (sigmaU):", value = 0.0014))
      ),
      fluidRow(
        column(6, numericInput("YR", "Restricted estimate (YR):", value = 0.0026)),
        column(6, numericInput("sigmaR", "Standard error (sigmaR):", value = 0.0009))
      ),
      numericInput("VUR", "Cov(YU,YR) (default assumes YR is efficient, and sets VUR to sigmaR^2):", value = 0.0009^2, step = 0.1),
      actionButton("compute", "Adapt")
    ),
    
    mainPanel(
      textOutput("corr_output"),
      textOutput("st_output"),
      textOutput("estimate_output")
      
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
  # Precompute the function
  st_thresholds_grid <- c(1.57482766198508, 1.55125499909256, 1.52742074353448, 1.50334882329505, 1.47906233490077, 
                          1.45458150029266, 1.4299270409132, 1.40512383967876, 1.38019303478176, 1.35515605260822, 
                          1.33003227704614, 1.30483757026634, 1.27959849366401, 1.25433502048039, 1.22906433002374, 
                          1.20379862855458, 1.17856184780166, 1.15337076343841, 1.12824236155291, 1.10319422804627, 
                          1.07824432631491, 1.05340561586754, 1.02870454439195, 1.00416561484403, 0.979813578603983, 
                          0.955673742841495, 0.931764912197942, 0.908115155200027, 0.884752852219577, 0.861705205296034, 
                          0.838984986407601, 0.816633434009214, 0.794675112145209, 0.7731297065466, 0.752040905293611, 
                          0.731416663001888, 0.711309030584423, 0.69172503640683, 0.672713172851095, 0.654291308290986, 
                          0.636491137084625, 0.619350676770425, 0.602881936878391, 0.58712468096191, 0.57210936993747, 
                          0.557856650375231, 0.544387216517243, 0.531735712033272, 0.519924940086376, 0.508976432315778, 
                          0.498906428201768, 0.489732214564344, 0.481477367692955, 0.474157387541963, 0.467784346318755, 
                          0.462372207805892, 0.457930286791032, 0.454466861185033, 0.451997110894451, 0.450509551844624)
  
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
  output$corr_output <- renderText({
    paste("The correlation coefficient between YU and (YR-YU) is ", round(corr,3),
          ", which implies the efficient estimate is ",GMM,
          " and the efficiency of YU relative to the efficient estimate is ", round(1-corr^2,3),".")
  })
  output$st_output <- renderText({
    paste("Given the relative efficiency, the soft threshold that achieves optimal adaptation is ", round(st,3),".")
  })
  
  if (-st < tO & tO < st) {
    output$estimate_output <- renderText({
      paste("Since (YR-YU)/std(YR-YU) = ",tO,
            " does not exceed the threshold, adaptation maintains to be the efficient estimate ",
            adaptive_st,".")
      })
  }
  if (tO <=-st ) {
    output$estimate_output <- renderText({
      paste("Since (YR-YU)/std(YR-YU) = ",tO,
            " is negative and below threshold, adaptation translates YU toward the efficient estimate and obtains ",
            adaptive_st,".")
    })
  }
  if (tO >= st ) {
    output$estimate_output <- renderText({
      paste("Since (YR-YU)/std(YR-YU) = ",tO,
            " is positive and above threshold, adaptation translates YU toward the efficient estimate and obtains ",
            adaptive_st,".")
    })
  }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
