library(shiny)
library(ggplot2)

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

  # Third section (m) -------------
  p("Calculate $m$: number of individuals per cluster", style = "color:blue; font-size:25px"),

  p("Input Parameters", style = "font-size:20px"),

  # Inputs at the top
  fluidRow(
    column(3,
           textInput("power", "Power:", value = ""),
           textInput("beta1", "Effect for Y1:", value = ""),
           textInput("beta2", "Effect for Y2:", value = ""),
           textInput("alpha", "Type I error:", value = "")
    ),
    column(3,
           textInput("K", "K (# clusters in treatment group):", value = ""),
           textInput("varY1", "Total variance of Y1:", value = ""),
           textInput("varY2", "Total variance of Y2:", value = ""),
           textInput("rho1", "Inter-subject between-endpoint ICC:", value = "")
    ),
    column(3,
           textInput("r", "Treatment allocation ratio:", value = ""),
           textInput("rho01", "ICC for Y1:", value = ""),
           textInput("rho02", "ICC for Y2:", value = ""),
           textInput("rho2", "Intra-subject between-endpoint ICC:", value = "")
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

  # Section results 3 (m) ---------------
  # Reactive expression to run the functions and create the data frame
  results3 <- eventReactive(input$calcButton3, {
    # Get user inputs
    r_input <- as.numeric(input$r)
    power_input <- as.numeric(input$power)
    K_input <- as.numeric(input$K)
    alpha_input <- as.numeric(input$alpha)
    beta1_input <- as.numeric(input$beta1)
    beta2_input <- as.numeric(input$beta2)
    varY1_input <- as.numeric(input$varY1)
    varY2_input <- as.numeric(input$varY2)
    rho01_input <- as.numeric(input$rho01)
    rho02_input <- as.numeric(input$rho02)
    rho1_input <- as.numeric(input$rho1)
    rho2_input <- as.numeric(input$rho2)

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

  # Output the bargraph 3 for m -------------
  output$bargraph3 <- renderPlot({
    # Get the results
    df <- results3()

    # Check if the data frame is empty
    if (nrow(df) > 0) {
      # Create the bar graph using ggplot2
      ggplot(df, aes(x = reorder(Function, Fill), y = Value, fill = as.factor(Fill))) +
        geom_bar(stat = "identity") +
        ylab("m") +
        xlab("Design Method") +
        geom_text(aes(label = Value), vjust = -0.5, size = 4) +
        ggtitle("Figure 3. Calculations for m - number of individuals in each cluster") +
        theme(legend.position = "none")
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
