library(shiny)
library(ggplot2)
library(latex2exp)

# Load your package
# library(yourpackage)  # Uncomment and replace 'yourpackage' with the actual package name

# Define UI for application
ui <- fluidPage(
  tags$head(
    tags$link(rel="stylesheet",
              href="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.css",
              integrity="sha384-dbVIfZGuN1Yq7/1Ocstc1lUEm+AT+/rCkibIcC/OmWo5f0EA48Vf8CytHzGrSwbQ",
              crossorigin="anonymous"),
    HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.js" integrity="sha384-2BKqo+exmr9su6dir+qCw08N2ZKRucY4PrGQPPWU1A7FtlCGjmEGFqXCv5nyM5Ij" crossorigin="anonymous"></script>'),
    HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous"></script>'),
    HTML('
    <script>
      document.addEventListener("DOMContentLoaded", function(){
        renderMathInElement(document.body, {
          delimiters: [{left: "$", right: "$", display: false}]
        });
      })
    </script>')
  ),

  titlePanel("Calculate study design specifications for cluster-randomized trials \nwith co-primary outcomes using `crt2power` package"),

  p("This ShinyApp lets the user calculate the number of clusters in the treatment group ($K$), cluster size ($m$), or statistical power ($\\pi$) from the user's desired input parameters. Calculations are done using the R package `crt2power`."),

  shinyjs::useShinyjs(),
  shinyjs::extendShinyjs(text = "shinyjs.refresh_page = function() { location.reload(); }", functions = "refresh_page"),
  actionButton("refresh", "Refresh Application"),

  # UI 1 (K) -------------
  p("Calculate Power: probability of correctly rejecting the null hypothesis", style = "color:blue; font-size:25px"),

  p("Input Parameters", style = "font-size:20px"),

  # Inputs at the top
  fluidRow(
    column(3,
           textInput("K_power", "$K_1$: Number of clusters in treatment group", value = ""),
           textInput("beta1_power", "$\\beta_1$: Effect for Y1", value = ""),
           textInput("beta2_power", "$\\beta_2$: Effect for Y2", value = ""),
           textInput("alpha_power", "$\\alpha$: Type I error", value = "")
    ),
    column(3,
           textInput("m_power", "$m$: Number of subjects in each cluster", value = ""),
           textInput("varY1_power", "Var$(Y_1)$: Total variance of outcome 1", value = ""),
           textInput("varY2_power", "Var$(Y_2)$: Total variance of outcome 2", value = ""),
           textInput("rho1_power", "$\\rho_1^{(1,2)}$: Inter-subject between-endpoint ICC", value = "")
    ),
    column(3,
           textInput("r_power", "$r = K_2 / K_1$: Treatment allocation ratio", value = ""),
           textInput("rho01_power", "$\\rho_0^{(1)}$: ICC for outcome 1", value = ""),
           textInput("rho02_power", "$\\rho_0^{(2)}$: ICC for outcome 2", value = ""),
           textInput("rho2_power", "$\\rho_2^{(1,2)}$: Intra-subject between-endpoint ICC", value = "")
    )
  ),

  # Action button to trigger the calculation
  fluidRow(
    column(12,
           actionButton("calcButton1", "Calculate")
    )
  ),

  # Output: Bargraph
  fluidRow(
    column(12,
           plotOutput("bargraph1")
    )
  ),

  # UI 2 (K) -------------
  p("Calculate $K$: number of clusters in treatment group", style = "color:blue; font-size:25px"),

  p("Input Parameters", style = "font-size:20px"),

  # Inputs at the top
  fluidRow(
    column(3,
           textInput("power_K", "$\\pi$: Statistical power", value = ""),
           textInput("beta1_K", "$\\beta_1$: Effect for Y1", value = ""),
           textInput("beta2_K", "$\\beta_2$: Effect for Y2", value = ""),
           textInput("alpha_K", "$\\alpha$: Type I error", value = "")
    ),
    column(3,
           textInput("m_K", "$m$: Number of subjects in each cluster", value = ""),
           textInput("varY1_K", "Var$(Y_1)$: Total variance of outcome 1", value = ""),
           textInput("varY2_K", "Var$(Y_2)$: Total variance of outcome 2", value = ""),
           textInput("rho1_K", "$\\rho_1^{(1,2)}$: Inter-subject between-endpoint ICC", value = "")
    ),
    column(3,
           textInput("r_K", "$r = K_2 / K_1$: Treatment allocation ratio", value = ""),
           textInput("rho01_K", "$\\rho_0^{(1)}$: ICC for outcome 1", value = ""),
           textInput("rho02_K", "$\\rho_0^{(2)}$: ICC for outcome 2", value = ""),
           textInput("rho2_K", "$\\rho_2^{(1,2)}$: Intra-subject between-endpoint ICC", value = "")
    )
  ),

  # Action button to trigger the calculation
  fluidRow(
    column(12,
           actionButton("calcButton2", "Calculate")
    )
  ),

  # Output: Bargraph
  fluidRow(
    column(12,
           plotOutput("bargraph2")
    )
  ),

  # UI 3 (m) -------------
  p("Calculate $m$: number of individuals per cluster", style = "color:blue; font-size:25px"),

  p("Input Parameters", style = "font-size:20px"),

  # Inputs at the top
  fluidRow(
    column(3,
           textInput("power_m", "$\\pi$: Statistical power", value = ""),
           textInput("beta1_m", "$\\beta_1$: Effect for Y1", value = ""),
           textInput("beta2_m", "$\\beta_2$: Effect for Y2", value = ""),
           textInput("alpha_m", "$\\alpha$: Type I error", value = "")
    ),
    column(3,
           textInput("K_m", "$K_1$: Number of clusters in treatment group", value = ""),
           textInput("varY1_m", "Var$(Y_1)$: Total variance of outcome 1", value = ""),
           textInput("varY2_m", "Var$(Y_2)$: Total variance of outcome 2", value = ""),
           textInput("rho1_m", "$\\rho_1^{(1,2)}$: Inter-subject between-endpoint ICC", value = "")
    ),
    column(3,
           textInput("r_m", "$r = K_2 / K_1$: Treatment allocation ratio", value = ""),
           textInput("rho01_m", "$\\rho_0^{(1)}$: ICC for outcome 1", value = ""),
           textInput("rho02_m", "$\\rho_0^{(2)}$: ICC for outcome 2", value = ""),
           textInput("rho2_m", "$\\rho_2^{(1,2)}$: Intra-subject between-endpoint ICC", value = "")
    )
  ),

  # Action button to trigger the calculation
  fluidRow(
    column(12,
           actionButton("calcButton3", "Calculate")
    )
  ),

  # Output: Bargraph
  fluidRow(
    column(12,
           plotOutput("bargraph3")
    )
  )
)

# Define server logic
server <- function(input, output) {

  observeEvent(input$refresh, {
    shinyjs::js$refresh_page()
  })

# Server 1 (Power) ---------------
  # Reactive expression to run the functions and create the data frame
  results1 <- eventReactive(input$calcButton1, {
    # Get user inputs
    r_input <- as.numeric(input$r_power)
    m_input <- as.numeric(input$m_power)
    K_input <- as.numeric(input$K_power)
    alpha_input <- as.numeric(input$alpha_power)
    beta1_input <- as.numeric(input$beta1_power)
    beta2_input <- as.numeric(input$beta2_power)
    varY1_input <- as.numeric(input$varY1_power)
    varY2_input <- as.numeric(input$varY2_power)
    rho01_input <- as.numeric(input$rho01_power)
    rho02_input <- as.numeric(input$rho02_power)
    rho1_input <- as.numeric(input$rho1_power)
    rho2_input <- as.numeric(input$rho2_power)

    r_input <- 1
    m_input <- 300
    K_input <- 15
    alpha_input <- 0.05
    beta1_input <- 0.1
    beta2_input <- 0.1
    varY1_input <- 0.23
    varY2_input <- 0.25
    rho01_input <- 0.025
    rho02_input <- 0.025
    rho1_input <- 0.01
    rho2_input <- 0.05

    # Ensure inputs are not empty
    if (!is.na(r_input) && !is.na(m_input) && !is.na(K_input) && !is.na(alpha_input) && !is.na(beta1_input) && !is.na(beta2_input) && !is.na(varY1_input) && !is.na(varY2_input) && !is.na(rho01_input) && !is.na(rho02_input) && !is.na(rho1_input) && !is.na(rho2_input)) {
      power1_bonf <- calc_pwr_pval_adj(K = K_input, m = m_input,
                                     alpha = alpha_input,
                                     beta1 = beta1_input, beta2 = beta2_input,
                                     varY1 = varY1_input, varY2 = varY2_input,
                                     rho01 = rho01_input, rho02 = rho02_input,
                                     rho2  = rho2_input, r = r_input)$`Final Power`[1]

      power1_sidak <- calc_pwr_pval_adj(K = K_input, m = m_input,
                                      alpha = alpha_input,
                                      beta1 = beta1_input, beta2 = beta2_input,
                                      varY1 = varY1_input, varY2 = varY2_input,
                                      rho01 = rho01_input, rho02 = rho02_input,
                                      rho2  = rho2_input, r = r_input)$`Final Power`[2]

      power1_dap <- calc_pwr_pval_adj(K = K_input, m = m_input,
                                alpha = alpha_input,
                                beta1 = beta1_input, beta2 = beta2_input,
                                varY1 = varY1_input, varY2 = varY2_input,
                                rho01 = rho01_input, rho02 = rho02_input,
                                rho2  = rho2_input, r = r_input)$`Final Power`[3]

      power2 <- calc_pwr_comb_outcome(K = K_input, m = m_input,
                                alpha = alpha_input, r = r_input,
                                beta1 = beta1_input, beta2 = beta2_input,
                                varY1 = varY1_input, varY2 = varY2_input,
                                rho01 = rho01_input, rho02 = rho02_input,
                                rho1 = rho1_input, rho2  = rho2_input)

      power3 <- calc_pwr_single_1dftest(m = m_input, K = K_input,
                                  alpha = alpha_input, r = r_input,
                                  beta1 = beta1_input, beta2 = beta2_input,
                                  varY1 = varY1_input, varY2 = varY2_input,
                                  rho01 = rho01_input, rho02 = rho02_input,
                                  rho1 = rho1_input, rho2  = rho2_input)

      power4_F <- calc_pwr_disj_2dftest(m = m_input, K = K_input,
                                  alpha = alpha_input, r = r_input,
                                  beta1 = beta1_input, beta2 = beta2_input,
                                  varY1 = varY1_input, varY2 = varY2_input,
                                  rho01 = rho01_input, rho02 = rho02_input,
                                  rho1 = rho1_input, rho2  = rho2_input,
                                  dist = "F")

      power4_Chi2 <- calc_pwr_disj_2dftest(m = m_input, K = K_input,
                                     alpha = alpha_input, r = r_input,
                                     beta1 = beta1_input, beta2 = beta2_input,
                                     varY1 = varY1_input, varY2 = varY2_input,
                                     rho01 = rho01_input, rho02 = rho02_input,
                                     rho1 = rho1_input, rho2  = rho2_input,
                                     dist = "Chi2")

      power5_T <- calc_pwr_conj_test(m = m_input, K = K_input,
                               alpha = alpha_input, r = r_input,
                               beta1 = beta1_input, beta2 = beta2_input,
                               varY1 = varY1_input, varY2 = varY2_input,
                               rho01 = rho01_input, rho02 = rho02_input,
                               rho1 = rho1_input, rho2  = rho2_input,
                               dist = "T")

      power5_MVN <- calc_pwr_conj_test(m = m_input, K = K_input,
                                 alpha = alpha_input, r = r_input,
                                 beta1 = beta1_input, beta2 = beta2_input,
                                 varY1 = varY1_input, varY2 = varY2_input,
                                 rho01 = rho01_input, rho02 = rho02_input,
                                 rho1 = rho1_input, rho2  = rho2_input,
                                 dist = "MVN")

      # Create a data frame with the results
      data.frame(
        Function = c("P-Value Adj. (Bonferonni)", "P-Value Adj. (Sidak)",
                     "P-Value Adj. (D/AP)", "Combined Outcomes",
                     "Single Weighted", "Disjunctive F-dist",
                     "Disjunctive Chi2", "Conjunctive T", "Conjunctive MVN"),
        Value = c(power1_bonf, power1_sidak, power1_dap,
                  power2, power3, power4_F, power4_Chi2,
                  power5_T, power5_MVN),
        Fill = c(1, 1, 1, 2, 3, 4, 4, 5, 5)
      )
    } else {
      data.frame(
        Function = character(0),
        Value = numeric(0)
      )
    }
  })

  # Output 1 (m) -------------
  output$bargraph1 <- renderPlot({
    # Get the results
    df <- results1()

    # Check if the data frame is empty
    if (nrow(df) > 0) {
      # Create the bar graph using ggplot2
      ggplot(df, aes(x = reorder(Function, Fill), y = Value)) +
        geom_bar(stat = "identity", fill = 'purple') +
        ylab("Power") +
        xlab("Design Method") +
        geom_text(aes(label = Value), vjust = -0.5, size = 4) +
        ggtitle("Figure 1. Results for statistical power")
    }
  })

# Server 2 (K) ---------------
  # Reactive expression to run the functions and create the data frame
  results2 <- eventReactive(input$calcButton2, {
    # Get user inputs
    r_input <- as.numeric(input$r_K)
    power_input <- as.numeric(input$power_K)
    m_input <- as.numeric(input$m_K)
    alpha_input <- as.numeric(input$alpha_K)
    beta1_input <- as.numeric(input$beta1_K)
    beta2_input <- as.numeric(input$beta2_K)
    varY1_input <- as.numeric(input$varY1_K)
    varY2_input <- as.numeric(input$varY2_K)
    rho01_input <- as.numeric(input$rho01_K)
    rho02_input <- as.numeric(input$rho02_K)
    rho1_input <- as.numeric(input$rho1_K)
    rho2_input <- as.numeric(input$rho2_K)

    # r_input <- 0.5
    # power_input <- 0.8
    # m_input <- 300
    # alpha_input <- 0.05
    # beta1_input <- 0.1
    # beta2_input <- 0.1
    # varY1_input <- 0.23
    # varY2_input <- 0.25
    # rho01_input <- 0.025
    # rho02_input <- 0.025
    # rho1_input <- 0.01
    # rho2_input <- 0.05

    # Ensure inputs are not empty
    if (!is.na(r_input) && !is.na(power_input) && !is.na(m_input) && !is.na(alpha_input) && !is.na(beta1_input) && !is.na(beta2_input) && !is.na(varY1_input) && !is.na(varY2_input) && !is.na(rho01_input) && !is.na(rho02_input) && !is.na(rho1_input) && !is.na(rho2_input)) {
      if(r_input == 1){
        K1_bonf_trt <- calc_K_pval_adj(m = m_input, power = power_input,
                                       alpha = alpha_input,
                                       beta1 = beta1_input, beta2 = beta2_input,
                                       varY1 = varY1_input, varY2 = varY2_input,
                                       rho01 = rho01_input, rho02 = rho02_input,
                                       rho2  = rho2_input,
                                       r = r_input)$`Final Treatment (K)`[[1]]
        K1_bonf_ctl <- calc_K_pval_adj(m = m_input, power = power_input,
                                       alpha = alpha_input,
                                       beta1 = beta1_input, beta2 = beta2_input,
                                       varY1 = varY1_input, varY2 = varY2_input,
                                       rho01 = rho01_input, rho02 = rho02_input,
                                       rho2  = rho2_input,
                                       r = r_input)$`Final Control (K)`[[1]]

        K1_sidak_trt <- calc_K_pval_adj(m = m_input, power = power_input,
                                        alpha = alpha_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho2  = rho2_input,
                                        r = r_input)$`Final Treatment (K)`[[2]]
        K1_sidak_ctl <- calc_K_pval_adj(m = m_input, power = power_input,
                                        alpha = alpha_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho2  = rho2_input,
                                        r = r_input)$`Final Control (K)`[[2]]

        K1_dap_trt <- calc_K_pval_adj(m = m_input, power = power_input,
                                      alpha = alpha_input,
                                      beta1 = beta1_input, beta2 = beta2_input,
                                      varY1 = varY1_input, varY2 = varY2_input,
                                      rho01 = rho01_input, rho02 = rho02_input,
                                      rho2  = rho2_input,
                                      r = r_input)$`Final Treatment (K)`[[3]]
        K1_dap_ctl <- calc_K_pval_adj(m = m_input, power = power_input,
                                      alpha = alpha_input,
                                      beta1 = beta1_input, beta2 = beta2_input,
                                      varY1 = varY1_input, varY2 = varY2_input,
                                      rho01 = rho01_input, rho02 = rho02_input,
                                      rho2  = rho2_input,
                                      r = r_input)$`Final Control (K)`[[3]]

        K2_trt <- calc_K_comb_outcome(m = m_input, power = power_input,
                                      alpha = alpha_input, r = r_input,
                                      beta1 = beta1_input, beta2 = beta2_input,
                                      varY1 = varY1_input, varY2 = varY2_input,
                                      rho01 = rho01_input, rho02 = rho02_input,
                                      rho1 = rho1_input,
                                      rho2 = rho2_input)$`Treatment (K)`[[1]]
        K2_ctl <- calc_K_comb_outcome(m = m_input, power = power_input,
                                      alpha = alpha_input, r = r_input,
                                      beta1 = beta1_input, beta2 = beta2_input,
                                      varY1 = varY1_input, varY2 = varY2_input,
                                      rho01 = rho01_input, rho02 = rho02_input,
                                      rho1 = rho1_input,
                                      rho2 = rho2_input)$`Control (K)`[[1]]

        K3_trt <- calc_K_single_1dftest(power = power_input, m = m_input,
                                        alpha = alpha_input, r = r_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho1 = rho1_input,
                                        rho2 = rho2_input)$`Treatment (K)`[[1]]
        K3_ctl <- calc_K_single_1dftest(power = power_input, m = m_input,
                                        alpha = alpha_input, r = r_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho1 = rho1_input,
                                        rho2 = rho2_input)$`Control (K)`[[1]]

        K4_F_trt <- calc_K_disj_2dftest(power = power_input, m = m_input,
                                        alpha = alpha_input, r = r_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho1 = rho1_input, rho2  = rho2_input,
                                        dist = "F")$`Treatment (K)`[[1]]
        K4_F_ctl <- calc_K_disj_2dftest(power = power_input, m = m_input,
                                        alpha = alpha_input, r = r_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho1 = rho1_input, rho2  = rho2_input,
                                        dist = "F")$`Control (K)`[[1]]

        K4_Chi2_trt <- calc_K_disj_2dftest(power = power_input, m = m_input,
                                           alpha = alpha_input, r = r_input,
                                           beta1 = beta1_input, beta2 = beta2_input,
                                           varY1 = varY1_input, varY2 = varY2_input,
                                           rho01 = rho01_input, rho02 = rho02_input,
                                           rho1 = rho1_input, rho2  = rho2_input,
                                           dist = "Chi2")$`Treatment (K)`[[1]]
        K4_Chi2_ctl <- calc_K_disj_2dftest(power = power_input, m = m_input,
                                           alpha = alpha_input, r = r_input,
                                           beta1 = beta1_input, beta2 = beta2_input,
                                           varY1 = varY1_input, varY2 = varY2_input,
                                           rho01 = rho01_input, rho02 = rho02_input,
                                           rho1 = rho1_input, rho2  = rho2_input,
                                           dist = "Chi2")$`Control (K)`[[1]]

        K5_T_trt <- calc_K_conj_test(power = power_input, m = m_input,
                                     alpha = alpha_input, r = r_input,
                                     beta1 = beta1_input, beta2 = beta2_input,
                                     varY1 = varY1_input, varY2 = varY2_input,
                                     rho01 = rho01_input, rho02 = rho02_input,
                                     rho1 = rho1_input, rho2  = rho2_input,
                                     dist = "T")$`Treatment (K)`[[1]]
        K5_T_ctl <- calc_K_conj_test(power = power_input, m = m_input,
                                     alpha = alpha_input, r = r_input,
                                     beta1 = beta1_input, beta2 = beta2_input,
                                     varY1 = varY1_input, varY2 = varY2_input,
                                     rho01 = rho01_input, rho02 = rho02_input,
                                     rho1 = rho1_input, rho2  = rho2_input,
                                     dist = "T")$`Control (K)`[[1]]

        K5_MVN_trt <- calc_K_conj_test(power = power_input, m = m_input,
                                       alpha = alpha_input, r = r_input,
                                       beta1 = beta1_input, beta2 = beta2_input,
                                       varY1 = varY1_input, varY2 = varY2_input,
                                       rho01 = rho01_input, rho02 = rho02_input,
                                       rho1 = rho1_input, rho2  = rho2_input,
                                       dist = "MVN")$`Treatment (K)`[[1]]
        K5_MVN_ctl <- calc_K_conj_test(power = power_input, m = m_input,
                                       alpha = alpha_input, r = r_input,
                                       beta1 = beta1_input, beta2 = beta2_input,
                                       varY1 = varY1_input, varY2 = varY2_input,
                                       rho01 = rho01_input, rho02 = rho02_input,
                                       rho1 = rho1_input, rho2  = rho2_input,
                                       dist = "MVN")$`Control (K)`[[1]]

      } else{
        K1_bonf_trt <- calc_K_pval_adj(m = m_input, power = power_input,
                                       alpha = alpha_input,
                                       beta1 = beta1_input, beta2 = beta2_input,
                                       varY1 = varY1_input, varY2 = varY2_input,
                                       rho01 = rho01_input, rho02 = rho02_input,
                                       rho2  = rho2_input,
                                       r = r_input)$`Final Treatment (K1)`[[1]]
        K1_bonf_ctl <- calc_K_pval_adj(m = m_input, power = power_input,
                                       alpha = alpha_input,
                                       beta1 = beta1_input, beta2 = beta2_input,
                                       varY1 = varY1_input, varY2 = varY2_input,
                                       rho01 = rho01_input, rho02 = rho02_input,
                                       rho2  = rho2_input,
                                       r = r_input)$`Final Control (K2)`[[1]]

        K1_sidak_trt <- calc_K_pval_adj(m = m_input, power = power_input,
                                        alpha = alpha_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho2  = rho2_input,
                                        r = r_input)$`Final Treatment (K1)`[[2]]
        K1_sidak_ctl <- calc_K_pval_adj(m = m_input, power = power_input,
                                        alpha = alpha_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho2  = rho2_input,
                                        r = r_input)$`Final Control (K2)`[[2]]

        K1_dap_trt <- calc_K_pval_adj(m = m_input, power = power_input,
                                      alpha = alpha_input,
                                      beta1 = beta1_input, beta2 = beta2_input,
                                      varY1 = varY1_input, varY2 = varY2_input,
                                      rho01 = rho01_input, rho02 = rho02_input,
                                      rho2  = rho2_input,
                                      r = r_input)$`Final Treatment (K1)`[[3]]
        K1_dap_ctl <- calc_K_pval_adj(m = m_input, power = power_input,
                                      alpha = alpha_input,
                                      beta1 = beta1_input, beta2 = beta2_input,
                                      varY1 = varY1_input, varY2 = varY2_input,
                                      rho01 = rho01_input, rho02 = rho02_input,
                                      rho2  = rho2_input,
                                      r = r_input)$`Final Control (K2)`[[3]]

        K2_trt <- calc_K_comb_outcome(m = m_input, power = power_input,
                                      alpha = alpha_input, r = r_input,
                                      beta1 = beta1_input, beta2 = beta2_input,
                                      varY1 = varY1_input, varY2 = varY2_input,
                                      rho01 = rho01_input, rho02 = rho02_input,
                                      rho1 = rho1_input,
                                      rho2 = rho2_input)$`Treatment (K1)`[[1]]
        K2_ctl <- calc_K_comb_outcome(m = m_input, power = power_input,
                                      alpha = alpha_input, r = r_input,
                                      beta1 = beta1_input, beta2 = beta2_input,
                                      varY1 = varY1_input, varY2 = varY2_input,
                                      rho01 = rho01_input, rho02 = rho02_input,
                                      rho1 = rho1_input,
                                      rho2 = rho2_input)$`Control (K2)`[[1]]

        K3_trt <- calc_K_single_1dftest(power = power_input, m = m_input,
                                        alpha = alpha_input, r = r_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho1 = rho1_input,
                                        rho2 = rho2_input)$`Treatment (K1)`[[1]]
        K3_ctl <- calc_K_single_1dftest(power = power_input, m = m_input,
                                        alpha = alpha_input, r = r_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho1 = rho1_input,
                                        rho2 = rho2_input)$`Control (K2)`[[1]]

        K4_F_trt <- calc_K_disj_2dftest(power = power_input, m = m_input,
                                        alpha = alpha_input, r = r_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho1 = rho1_input, rho2  = rho2_input,
                                        dist = "F")$`Treatment (K1)`[[1]]
        K4_F_ctl <- calc_K_disj_2dftest(power = power_input, m = m_input,
                                        alpha = alpha_input, r = r_input,
                                        beta1 = beta1_input, beta2 = beta2_input,
                                        varY1 = varY1_input, varY2 = varY2_input,
                                        rho01 = rho01_input, rho02 = rho02_input,
                                        rho1 = rho1_input, rho2  = rho2_input,
                                        dist = "F")$`Control (K2)`[[1]]

        K4_Chi2_trt <- calc_K_disj_2dftest(power = power_input, m = m_input,
                                           alpha = alpha_input, r = r_input,
                                           beta1 = beta1_input, beta2 = beta2_input,
                                           varY1 = varY1_input, varY2 = varY2_input,
                                           rho01 = rho01_input, rho02 = rho02_input,
                                           rho1 = rho1_input, rho2  = rho2_input,
                                           dist = "Chi2")$`Treatment (K1)`[[1]]
        K4_Chi2_ctl <- calc_K_disj_2dftest(power = power_input, m = m_input,
                                           alpha = alpha_input, r = r_input,
                                           beta1 = beta1_input, beta2 = beta2_input,
                                           varY1 = varY1_input, varY2 = varY2_input,
                                           rho01 = rho01_input, rho02 = rho02_input,
                                           rho1 = rho1_input, rho2  = rho2_input,
                                           dist = "Chi2")$`Control (K2)`[[1]]

        K5_T_trt <- calc_K_conj_test(power = power_input, m = m_input,
                                     alpha = alpha_input, r = r_input,
                                     beta1 = beta1_input, beta2 = beta2_input,
                                     varY1 = varY1_input, varY2 = varY2_input,
                                     rho01 = rho01_input, rho02 = rho02_input,
                                     rho1 = rho1_input, rho2  = rho2_input,
                                     dist = "T")$`Treatment (K1)`[[1]]
        K5_T_ctl <- calc_K_conj_test(power = power_input, m = m_input,
                                     alpha = alpha_input, r = r_input,
                                     beta1 = beta1_input, beta2 = beta2_input,
                                     varY1 = varY1_input, varY2 = varY2_input,
                                     rho01 = rho01_input, rho02 = rho02_input,
                                     rho1 = rho1_input, rho2  = rho2_input,
                                     dist = "T")$`Control (K2)`[[1]]

        K5_MVN_trt <- calc_K_conj_test(power = power_input, m = m_input,
                                       alpha = alpha_input, r = r_input,
                                       beta1 = beta1_input, beta2 = beta2_input,
                                       varY1 = varY1_input, varY2 = varY2_input,
                                       rho01 = rho01_input, rho02 = rho02_input,
                                       rho1 = rho1_input, rho2  = rho2_input,
                                       dist = "MVN")$`Treatment (K1)`[[1]]
        K5_MVN_ctl <- calc_K_conj_test(power = power_input, m = m_input,
                                       alpha = alpha_input, r = r_input,
                                       beta1 = beta1_input, beta2 = beta2_input,
                                       varY1 = varY1_input, varY2 = varY2_input,
                                       rho01 = rho01_input, rho02 = rho02_input,
                                       rho1 = rho1_input, rho2  = rho2_input,
                                       dist = "MVN")$`Control (K2)`[[1]]
      }

      # Create a data frame with the results
      data.frame(
        Function = c("P-Value Adj. (Bonferonni)", "P-Value Adj. (Sidak)",
                     "P-Value Adj. (D/AP)", "Combined Outcomes",
                     "Single Weighted", "Disjunctive F-dist",
                     "Disjunctive Chi2", "Conjunctive T", "Conjunctive MVN",
                     "P-Value Adj. (Bonferonni)", "P-Value Adj. (Sidak)",
                     "P-Value Adj. (D/AP)", "Combined Outcomes",
                     "Single Weighted", "Disjunctive F-dist",
                     "Disjunctive Chi2", "Conjunctive T", "Conjunctive MVN"),
        Value = c(K1_bonf_trt, K1_sidak_trt, K1_dap_trt,
                  K2_trt, K3_trt,
                  K4_F_trt, K4_Chi2_trt,
                  K5_T_trt, K5_MVN_trt,
                  K1_bonf_ctl, K1_sidak_ctl, K1_dap_ctl,
                  K2_ctl, K3_ctl,
                  K4_F_ctl, K4_Chi2_ctl,
                  K5_T_ctl, K5_MVN_ctl),
        Group = c("Treatment", "Treatment", "Treatment",
                  "Treatment", "Treatment", "Treatment",
                  "Treatment", "Treatment", "Treatment",
                  "Control", "Control", "Control", "Control", "Control",
                  "Control", "Control", "Control", "Control"),
        Order = c(1, 1, 1, 2, 3, 4, 4, 5, 5,
                  1, 1, 1, 2, 3, 4, 4, 5, 5)
      )
    } else {
      data.frame(
        Function = character(0),
        Value = numeric(0)
      )
    }
  })

# Output 2 (K) -------------
  output$bargraph2 <- renderPlot({
    # Get the results
    df <- results2()

    # Check if the data frame is empty
    if (nrow(df) > 0) {
      # Create the bar graph using ggplot2
      print(ggplot(df, aes(x = reorder(Function, Order),
                           y = Value, fill = Group)) +
              geom_bar(stat = "identity", position = 'dodge') +
              ylab("Number of clusters") +
              xlab("Design Method") +
              geom_text(aes(label = Value),
                        position = position_dodge(width = 0.9),
                        vjust = -0.5, size = 4) +
              ggtitle("Figure 2. Results for number of clusters in treatment group (K1) and control group (K2)"))
    }
  })

# Server 3 (m) ---------------
  # Reactive expression to run the functions and create the data frame
  results3 <- eventReactive(input$calcButton3, {
    # Get user inputs
    r_input <- as.numeric(input$r_m)
    power_input <- as.numeric(input$power_m)
    K_input <- as.numeric(input$K_m)
    alpha_input <- as.numeric(input$alpha_m)
    beta1_input <- as.numeric(input$beta1_m)
    beta2_input <- as.numeric(input$beta2_m)
    varY1_input <- as.numeric(input$varY1_m)
    varY2_input <- as.numeric(input$varY2_m)
    rho01_input <- as.numeric(input$rho01_m)
    rho02_input <- as.numeric(input$rho02_m)
    rho1_input <- as.numeric(input$rho1_m)
    rho2_input <- as.numeric(input$rho2_m)

    # r_input <- 0.5
    # power_input <- 0.8
    # K_input <- 30
    # alpha_input <- 0.05
    # beta1_input <- 0.1
    # beta2_input <- 0.1
    # varY1_input <- 0.23
    # varY2_input <- 0.25
    # rho01_input <- 0.025
    # rho02_input <- 0.025
    # rho1_input <- 0.01
    # rho2_input <- 0.05

    # Ensure inputs are not empty
    if (!is.na(r_input) && !is.na(power_input) && !is.na(K_input) && !is.na(alpha_input) && !is.na(beta1_input) && !is.na(beta2_input) && !is.na(varY1_input) && !is.na(varY2_input) && !is.na(rho01_input) && !is.na(rho02_input) && !is.na(rho1_input) && !is.na(rho2_input)) {
      m1_bonf <- calc_m_pval_adj(K = K_input, power = power_input,
                                 alpha = alpha_input,
                                 beta1 = beta1_input, beta2 = beta2_input,
                                 varY1 = varY1_input, varY2 = varY2_input,
                                 rho01 = rho01_input, rho02 = rho02_input,
                                 rho2  = rho2_input, r = r_input)$`Final m`[1]

      m1_sidak <- calc_m_pval_adj(K = K_input, power = power_input,
                                  alpha = alpha_input,
                                  beta1 = beta1_input, beta2 = beta2_input,
                                  varY1 = varY1_input, varY2 = varY2_input,
                                  rho01 = rho01_input, rho02 = rho02_input,
                                  rho2  = rho2_input, r = r_input)$`Final m`[2]

      m1_dap <- calc_m_pval_adj(K = K_input, power = power_input,
                                alpha = alpha_input,
                                beta1 = beta1_input, beta2 = beta2_input,
                                varY1 = varY1_input, varY2 = varY2_input,
                                rho01 = rho01_input, rho02 = rho02_input,
                                rho2  = rho2_input, r = r_input)$`Final m`[3]

      m2 <- calc_m_comb_outcome(K = K_input, power = power_input,
                                alpha = alpha_input, r = r_input,
                                beta1 = beta1_input, beta2 = beta2_input,
                                varY1 = varY1_input, varY2 = varY2_input,
                                rho01 = rho01_input, rho02 = rho02_input,
                                rho1 = rho1_input, rho2  = rho2_input)

      m3 <- calc_m_single_1dftest(power = power_input, K = K_input,
                                  alpha = alpha_input, r = r_input,
                                  beta1 = beta1_input, beta2 = beta2_input,
                                  varY1 = varY1_input, varY2 = varY2_input,
                                  rho01 = rho01_input, rho02 = rho02_input,
                                  rho1 = rho1_input, rho2  = rho2_input)

      m4_F <- calc_m_disj_2dftest(power = power_input, K = K_input,
                                  alpha = alpha_input, r = r_input,
                                  beta1 = beta1_input, beta2 = beta2_input,
                                  varY1 = varY1_input, varY2 = varY2_input,
                                  rho01 = rho01_input, rho02 = rho02_input,
                                  rho1 = rho1_input, rho2  = rho2_input,
                                  dist = "F")

      m4_Chi2 <- calc_m_disj_2dftest(power = power_input, K = K_input,
                                     alpha = alpha_input, r = r_input,
                                     beta1 = beta1_input, beta2 = beta2_input,
                                     varY1 = varY1_input, varY2 = varY2_input,
                                     rho01 = rho01_input, rho02 = rho02_input,
                                     rho1 = rho1_input, rho2  = rho2_input,
                                     dist = "Chi2")

      m5_T <- calc_m_conj_test(power = power_input, K = K_input,
                               alpha = alpha_input, r = r_input,
                               beta1 = beta1_input, beta2 = beta2_input,
                               varY1 = varY1_input, varY2 = varY2_input,
                               rho01 = rho01_input, rho02 = rho02_input,
                               rho1 = rho1_input, rho2  = rho2_input,
                               dist = "T")

      m5_MVN <- calc_m_conj_test(power = power_input, K = K_input,
                                 alpha = alpha_input, r = r_input,
                                 beta1 = beta1_input, beta2 = beta2_input,
                                 varY1 = varY1_input, varY2 = varY2_input,
                                 rho01 = rho01_input, rho02 = rho02_input,
                                 rho1 = rho1_input, rho2  = rho2_input,
                                 dist = "MVN")

      # Create a data frame with the results
      data.frame(
        Function = c("P-Value Adj. (Bonferonni)", "P-Value Adj. (Sidak)",
                     "P-Value Adj. (D/AP)", "Combined Outcomes",
                     "Single Weighted", "Disjunctive F-dist",
                     "Disjunctive Chi2", "Conjunctive T", "Conjunctive MVN"),
        Value = c(m1_bonf, m1_sidak, m1_dap, m2, m3, m4_F, m4_Chi2, m5_T, m5_MVN),
        Fill = c(1, 1, 1, 2, 3, 4, 4, 5, 5)
      )
    } else {
      data.frame(
        Function = character(0),
        Value = numeric(0)
      )
    }
  })

# Output 3 (m) -------------
  output$bargraph3 <- renderPlot({
    # Get the results
    df <- results3()

    # Check if the data frame is empty
    if (nrow(df) > 0) {
      # Create the bar graph using ggplot2
      ggplot(df, aes(x = reorder(Function, Fill), y = Value)) +
        geom_bar(stat = "identity", fill = 'purple') +
        ylab("m") +
        xlab("Design Method") +
        geom_text(aes(label = Value), vjust = -0.5, size = 4) +
        ggtitle("Figure 3. Results for number of individuals per cluster (m)")
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
