library(shiny)
library(ggplot2)
library(latex2exp)
library(bslib)
library(shinydashboard)
library(crt2power)
library(tibble)
library(DT)
require(tidyverse)

# Load your package
# library(yourpackage)  # Uncomment and replace 'yourpackage' with the actual package name

# Define UI for application
ui <- page_navbar(
  tags$head(
    tags$link(rel="stylesheet",
              href="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.css",
              integrity="sha384-dbVIfZGuN1Yq7/1Ocstc1lUEm+AT+/rCkibIcC/OmWo5f0EA48Vf8CytHzGrSwbQ",
              crossorigin="anonymous"),
    HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.js" integrity="sha384-2BKqo+exmr9su6dir+qCw08N2ZKRucY4PrGQPPWU1A7FtlCGjmEGFqXCv5nyM5Ij" crossorigin="anonymous"></script>'),
    HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous"></script>'),
    HTML('<script> document.addEventListener("DOMContentLoaded", function(){
          renderMathInElement(document.body, {
          delimiters: [{left: "$", right: "$", display: false}]
        });
      })
    </script>'),
    #tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$style(HTML("
      .my_table .table > thead > tr > th {
        background-color: #6495ED;
        color: white;
        padding: 8px; /* Adjust padding as needed */
        text-align: left;
        vertical-align: middle;
        border-top: 2px solid #6495ED;
        border-left: 2px solid #6495ED; /* Left border for header */;
        border-right: 2px solid #6495ED; /* Right border for header */
        border-bottom: 2px solid #6495ED; /* Right border for header */
      }

      .my_table .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
        padding: 4px;
        text-align: left;
        line-height: 1.42857143;
        vertical-align: middle;
        border-top: 2px solid #6495ED;
        border-left: 2px solid #6495ED;
        border-right: 2px solid #6495ED;
        border-bottom: 2px solid #6495ED;
      }
    ")),

    tags$style(HTML("
    .result_table {
        display: flex; /* Use flexbox to center the table */
        justify-content: center; /* Center horizontally */
        margin-top: 20px; /* Add some top margin for spacing */
      }

      .result_table .table > thead > tr > th {
        background-color: #6495ED;
        color: white;
        border-top: 2px solid #6495ED;
        border-left: 2px solid #6495ED; /* Left border for header */;
        border-right: 2px solid #6495ED; /* Right border for header */
        border-bottom: 2px solid #6495ED; /* Right border for header */
      }

      .result_table .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
        border-top: 2px solid #6495ED;
        border-left: 2px solid #6495ED;
        border-right: 2px solid #6495ED;
        border-bottom: 2px solid #6495ED;
      }
    "))
  ),

  title = "`crt2power` Application",
  bg = "#FF69B4",
  inverse = TRUE,

  # Overview -------------
  nav_panel(title = "Overview",
            titlePanel(h1("Power and sample size calculator for cluster-randomized trials with two co-primary outcomes", align = "center")),
            fluidRow(
              column(4,
                     card(
                       card_header("Overview"),
                       "This ShinyApp lets the user calculate the number of clusters in the treatment group, cluster size, or statistical power from the user's desired input parameters. Calculations are done using the R package `crt2power`."
                       ),
                     card(
                       card_header("Optimizing Your Experience"),
                       "For the best viewing experience when using this application on your web browser, we recommend setting the zoom to 80%. All calculations should take less than a minute. If you are waiting on a calculation, it may be because the input parameters you have provided result in very low power, or too high of sample size requirements. When this happens, please use the `Refresh Application` button to reload the application."
                     ),
                     card(
                       card_header("Contact"),
                       p("Please see the links tab at the upper right for more information.",
                         tags$br(),
                         tags$b("Author: "), "Melody Owen",
                         tags$br(),
                         tags$b("Email: "), "melody.owen@yale.edu",
                         tags$br(),
                         tags$b("Affiliation: "), "Yale University, Department of Biostatistics, Center for Methods in Implementation and Prevention Science",
                         tags$br(),
                         tags$b("Acknowledgements: "), "Thank you to my advisors - Dr. Donna Spiegelman, Dr. Fan Li, and Dr. Laura Forastiere, and thank you to Yale for supporting this research"
                       )
                       )
                     ),
              column(8,
                     card(
                       card_header("Guide to Input Parameters"),
                       withMathJax(),
                       tags$div(class = "my_table", tableOutput('overviewTable')),
                       p("1. This assumes equal treatment allocation. To be more precise, we often refer to the number of clusters in the treatment group as $K_1$, and the number of clusters in the control group as $K_2$.",
                         tags$br(),
                         "2. Note that not all design methods make use of every parameter listed above.")
                       )
                     )
              ),
            fluid = TRUE), # End nav_panel()

  # UI 1 (power) -------------
  nav_panel(title = "Calculate Power",
            titlePanel("Calculate Power ($\\pi$): probability of correctly rejecting the null hypothesis"),
            p("This ShinyApp lets the user calculate the number of clusters in the treatment group, cluster size, or statistical power from the user's desired input parameters. Calculations are done using the R package `crt2power`. Enter your desired parameters in the boxes below, and click `Calculate` to generate the results for all five study design methods."),

            fluidRow( # Main row to contain both input columns and the bar graph

              column(4, # Input parameters and calculate button
                     card(
                       card_header("Input Parameters"),
                       fluidRow( # Inputs at the top
                         column(6,
                                textInput("K_power", "$K_1$: Clusters in treatment group", value = ""),
                                textInput("r_power", "$r = K_2 / K_1$: Allocation ratio", value = ""),
                                textInput("beta1_power", "$\\beta_1$: Effect for Y1", value = ""),
                                textInput("varY1_power", "Var$(Y_1)$: Total variance of $Y_1$", value = ""),
                                textInput("rho01_power", "$\\rho_0^{(1)}$: ICC for $Y_1$", value = ""),
                                textInput("rho1_power", "$\\rho_1^{(1,2)}$: Inter-subject between-endpoint ICC", value = "")
                                ),
                         column(6,
                                textInput("m_power", "$m$: Subjects in each cluster", value = ""),
                                textInput("alpha_power", "$\\alpha$: Type I error", value = ""),
                                textInput("beta2_power", "$\\beta_2$: Effect for Y2", value = ""),
                                textInput("varY2_power", "Var$(Y_2)$: Total variance of $Y_2$", value = ""),
                                textInput("rho02_power", "$\\rho_0^{(2)}$: ICC for $Y_2$", value = ""),
                                textInput("rho2_power", "$\\rho_2^{(1,2)}$: Intra-subject between-endpoint ICC", value = "")
                                )
                         ),
                       # Action button to trigger the calculation
                       fluidRow(
                         column(3, actionButton("calcButton1", "Calculate")),
                         column(5,
                                shinyjs::useShinyjs(),
                                shinyjs::extendShinyjs(text = "shinyjs.refresh_page = function() { location.reload(); }", functions = "refresh_page"),
                                actionButton("refresh", "Refresh Application")
                                )
                         ),
                       fluidRow = TRUE
                       )
                     ),

              # Bar graph to the right
              column(8,

                     navset_card_underline(
                       title = "Visualizations",
                       # Panel with plot ----
                       nav_panel("Plot", plotOutput("bargraph1", height = "587px")),

                       # Panel with table ----
                       #nav_panel("Table", tableOutput("table1"))

                       nav_panel("Table", tags$div(
                         class="result_table",
                         tableOutput('table1'))

                       )

                     )

              )
            ),
            fluid = TRUE),

  # UI 2 (K) -------------
  nav_panel(title = "Calculate K",
            titlePanel("Calculate $K$: number of clusters in treatment group"),
            p("This ShinyApp lets the user calculate the number of clusters in the treatment group, cluster size, or statistical power from the user's desired input parameters. Calculations are done using the R package `crt2power`. Enter your desired parameters in the boxes below, and click `Calculate` to generate the results for all five study design methods."),

            # Main row to contain both input columns and the bar graph
            fluidRow(
              # Input parameters and calculate button
              column(4,
                     card(card_header("Input Parameters"),
                     # Inputs at the top
                     fluidRow(
                       column(6,
                              textInput("power_K", "$\\pi$: Statistical power", value = ""),
                              textInput("r_K", "$r = K_2 / K_1$: Allocation ratio", value = ""),
                              textInput("beta1_K", "$\\beta_1$: Effect for Y1", value = ""),
                              textInput("varY1_K", "Var$(Y_1)$: Total variance of $Y_1$", value = ""),
                              textInput("rho01_K", "$\\rho_0^{(1)}$: ICC for $Y_1$", value = ""),
                              textInput("rho1_K", "$\\rho_1^{(1,2)}$: Inter-subject between-endpoint ICC", value = "")
                       ),
                       column(6,
                              textInput("m_K", "$m$: Subjects in each cluster", value = ""),
                              textInput("alpha_K", "$\\alpha$: Type I error", value = ""),
                              textInput("beta2_K", "$\\beta_2$: Effect for Y2", value = ""),
                              textInput("varY2_K", "Var$(Y_2)$: Total variance of $Y_2$", value = ""),
                              textInput("rho02_K", "$\\rho_0^{(2)}$: ICC for $Y_2$", value = ""),
                              textInput("rho2_K", "$\\rho_2^{(1,2)}$: Intra-subject between-endpoint ICC", value = "")
                       )
                     ),
                     # Action button to trigger the calculation
                     fluidRow(
                       column(3,
                              actionButton("calcButton2", "Calculate")
                       ),
                       column(5,
                              shinyjs::useShinyjs(),
                              shinyjs::extendShinyjs(text = "shinyjs.refresh_page = function() { location.reload(); }", functions = "refresh_page"),
                              actionButton("refresh", "Refresh Application")
                       )
                     ),
                     fluidRow = TRUE
                     )
              ),
              # Bar graph to the right
              column(8,
                     navset_card_underline(
                       title = "Visualizations",
                       # Panel with plot ----
                       nav_panel("Plot", plotOutput("bargraph2", height = "587px")),

                       # Panel with table ----
                       #nav_panel("Table", tableOutput("table2"))

                       nav_panel("Table", tags$div(
                         class="result_table",
                         tableOutput('table2'))

                       )
                     )
              )
            ),
            fluid = TRUE),

  # UI 3 (m) -------------
  nav_panel(title = "Calculate m",
            titlePanel("Calculate $m$: number of individuals per cluster"),
            p("This ShinyApp lets the user calculate the number of clusters in the treatment group, cluster size, or statistical power from the user's desired input parameters. Calculations are done using the R package `crt2power`. Enter your desired parameters in the boxes below, and click `Calculate` to generate the results for all five study design methods."),

            # Main row to contain both input columns and the bar graph
            fluidRow(
              # Input parameters and calculate button
              column(4,
                     card(card_header("Input Parameters"),
                     # Inputs at the top
                     fluidRow(
                       column(6,
                              textInput("power_m", "$\\pi$: Statistical power", value = ""),
                              textInput("r_m", "$r = K_2 / K_1$: Allocation ratio", value = ""),
                              textInput("beta1_m", "$\\beta_1$: Effect for Y1", value = ""),
                              textInput("varY1_m", "Var$(Y_1)$: Total variance of $Y_1$", value = ""),
                              textInput("rho01_m", "$\\rho_0^{(1)}$: ICC for $Y_1$", value = ""),
                              textInput("rho1_m", "$\\rho_1^{(1,2)}$: Inter-subject between-endpoint ICC", value = "")
                       ),
                       column(6,
                              textInput("K_m", "$K_1$: Clusters in treatment group", value = ""),
                              textInput("alpha_m", "$\\alpha$: Type I error", value = ""),
                              textInput("beta2_m", "$\\beta_2$: Effect for Y2", value = ""),
                              textInput("varY2_m", "Var$(Y_2)$: Total variance of $Y_2$", value = ""),
                              textInput("rho02_m", "$\\rho_0^{(2)}$: ICC for $Y_2$", value = ""),
                              textInput("rho2_m", "$\\rho_2^{(1,2)}$: Intra-subject between-endpoint ICC", value = "")
                       )
                     ),
                     # Action button to trigger the calculation
                     fluidRow(
                       column(3,
                              actionButton("calcButton3", "Calculate")
                       ),
                       column(5,
                              shinyjs::useShinyjs(),
                              shinyjs::extendShinyjs(text = "shinyjs.refresh_page = function() { location.reload(); }", functions = "refresh_page"),
                              actionButton("refresh", "Refresh Application")
                       )
                     ),
                     fluidRow = TRUE
                     )
              ),
              # Bar graph to the right
              column(8,
                     navset_card_underline(
                       title = "Visualizations",
                       # Panel with plot ----
                       nav_panel("Plot", plotOutput("bargraph3", height = "587px")),

                       nav_panel("Table", tags$div(
                         class = "result_table",
                         tableOutput('table3'))

                       )
                     )
              )
            ),

            # shinyjs::useShinyjs(),
            # shinyjs::extendShinyjs(text = "shinyjs.refresh_page = function() { location.reload(); }", functions = "refresh_page"),
            # actionButton("refresh", "Refresh Application"),

            fluid = TRUE
            ), # End nav_panel()

  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("GitHub Repository", href = "https://github.com/melodyaowen/crt2power")),
    nav_item(tags$a("CRAN `crt2power` Manual", href = "https://shiny.posit.co"))
  ),


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

    # r_input <- 1
    # m_input <- 300
    # K_input <- 15
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
        Function = c("P-Value Adj. (Bonf.)", "P-Value Adj. (Sidak)",
                     "P-Value Adj. (D/AP)", "Combined Outcomes",
                     "Single Weighted", "Disjunctive F-dist",
                     "Disjunctive Chi2", "Conjunctive T-dist", "Conjunctive MVN"),
        TableLabel = c("P-Value Adjustment (Bonferroni)", "P-Value Adjustment (Sidak)",
                       "P-Value Adjustment (D/AP)", "Combined/Composite Outcomes Approach",
                       "Single 1-DF Weighted Test", "Disjunctive 2-DF Test (F Distribution)",
                       "Disjunctive 2-DF Test (Chi-2 Distribution)",
                       "Conjunctive IU Test (T Distribution)", "Conjunctive IU Test (Multivariate Normal Distribution)"),
        MethodIndex = c("Method 1", "Method 1", "Method 1",
                        "Method 2", "Method 3", "Method 4", "Method 4",
                        "Method 5", "Method 5"),
        Value = c(power1_bonf, power1_sidak, power1_dap,
                  power2, power3, power4_F, power4_Chi2,
                  power5_T, power5_MVN),
        Fill = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
      )
    } else {
      data.frame(
        Function = character(0),
        Value = numeric(0)
      )
    }
  })

  # Output 1 (power) -------------
  output$bargraph1 <- renderPlot({
    # Get the results
    df <- results1()

    # Check if the data frame is empty
    if (nrow(df) > 0) {
      # Create the bar graph using ggplot2
      ggplot(df, aes(x = reorder(Function, Fill), y = Value)) +
        geom_bar(stat = "identity", fill = '#6495ED') +
        ylab("Power") +
        xlab("Design Method") +
        ggtitle("Statistical Power Results") +
        geom_text(aes(label = Value), vjust = -0.5, size = 5) +
        theme(text = element_text(size = 20),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.line = element_line(colour = "black", size = 0.75,
                                       linetype = "solid")) +
        scale_y_continuous(limits = c(0, 1.1),
                           breaks = c(0, 0.25, 0.5, 0.75, 1))
    }
  })

  output$table1 <- renderTable({
    myTable1 <- results1() %>%
      dplyr::mutate(`Statistical Power` = paste0(round(Value*100, 2), "%")) %>%
      dplyr::select(`Method Group` = MethodIndex,
                    `Design Method` = TableLabel,
                    `Statistical Power`)


  }, sanitize.text.function = function(x) x)

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
        Function = c("P-Value Adj. (Bonf.)", "P-Value Adj. (Sidak)",
                     "P-Value Adj. (D/AP)", "Combined Outcomes",
                     "Single Weighted", "Disjunctive F-dist",
                     "Disjunctive Chi2", "Conjunctive T-dist", "Conjunctive MVN",
                     "P-Value Adj. (Bonf.)", "P-Value Adj. (Sidak)",
                     "P-Value Adj. (D/AP)", "Combined Outcomes",
                     "Single Weighted", "Disjunctive F-dist",
                     "Disjunctive Chi2", "Conjunctive T-dist", "Conjunctive MVN"),
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
        TableLabel = c("P-Value Adjustment (Bonferroni)", "P-Value Adjustment (Sidak)",
                       "P-Value Adjustment (D/AP)", "Combined/Composite Outcomes Approach",
                       "Single 1-DF Weighted Test", "Disjunctive 2-DF Test (F Distribution)",
                       "Disjunctive 2-DF Test (Chi-2 Distribution)",
                       "Conjunctive IU Test (T Distribution)",
                       "Conjunctive IU Test (Multivariate Normal Distribution)",
                       "P-Value Adjustment (Bonferroni)", "P-Value Adjustment (Sidak)",
                       "P-Value Adjustment (D/AP)", "Combined/Composite Outcomes Approach",
                       "Single 1-DF Weighted Test", "Disjunctive 2-DF Test (F Distribution)",
                       "Disjunctive 2-DF Test (Chi-2 Distribution)",
                       "Conjunctive IU Test (T Distribution)",
                       "Conjunctive IU Test (Multivariate Normal Distribution)"),
        MethodIndex = c("Method 1", "Method 1", "Method 1",
                        "Method 2", "Method 3", "Method 4", "Method 4",
                        "Method 5", "Method 5",
                        "Method 1", "Method 1", "Method 1",
                        "Method 2", "Method 3", "Method 4", "Method 4",
                        "Method 5", "Method 5"),
        Order = c(1, 2, 3, 4, 5, 6, 7, 8, 9,
                  1, 2, 3, 4, 5, 6, 7, 8, 9)
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
                        vjust = -0.5, size = 5) +
              ggtitle("Number of Clusters in Treatment Group (K1) and Control Group (K2) Results") +
              theme(text = element_text(size = 20),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    axis.line = element_line(colour = "black", size = 0.75,
                                             linetype = "solid"))
              )
    }
  })

  output$table2 <- renderTable({
    myTable3 <- results2() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Value = as.integer(Value)) %>%
      dplyr::arrange(Order) %>%
      tidyr::pivot_wider(names_from = Group, values_from = Value) %>%
      dplyr::select(`Method Group` = MethodIndex,
                    `Design Method` = TableLabel,
                    `Tretment (K1)` = Treatment,
                    `Control (K2)` = Control)
  }, sanitize.text.function = function(x) x)

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
        Function = c("P-Value Adj. (Bonf.)", "P-Value Adj. (Sidak)",
                     "P-Value Adj. (D/AP)", "Combined Outcomes",
                     "Single Weighted", "Disjunctive F-dist",
                     "Disjunctive Chi2", "Conjunctive T-dist", "Conjunctive MVN"),
        TableLabel = c("P-Value Adjustment (Bonferroni)", "P-Value Adjustment (Sidak)",
                       "P-Value Adjustment (D/AP)", "Combined/Composite Outcomes Approach",
                       "Single 1-DF Weighted Test", "Disjunctive 2-DF Test (F Distribution)",
                       "Disjunctive 2-DF Test (Chi-2 Distribution)",
                       "Conjunctive IU Test (T Distribution)", "Conjunctive IU Test (Multivariate Normal Distribution)"),
        MethodIndex = c("Method 1", "Method 1", "Method 1",
                        "Method 2", "Method 3", "Method 4", "Method 4",
                        "Method 5", "Method 5"),
        Value = c(m1_bonf, m1_sidak, m1_dap, m2, m3, m4_F, m4_Chi2, m5_T, m5_MVN),
        Fill = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
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
        geom_bar(stat = "identity", fill = '#6495ED') +
        ylab("Individuals per Cluster (m)") +
        xlab("Design Method") +
        geom_text(aes(label = Value), vjust = -0.5, size = 5) +
        ggtitle("Number of Individuals Per Cluster (m) Results") +
        theme(text = element_text(size = 20),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.line = element_line(colour = "black", size = 0.75,
                                       linetype = "solid"))

    }
  })

  output$table3 <- renderTable({
    myTable3 <- results3() %>%
      dplyr::mutate(Value = as.integer(Value)) %>%
      dplyr::select(`Method Group` = MethodIndex,
                    `Design Method` = TableLabel,
                    `Individuals Per Cluster (m)` = Value)

    }, sanitize.text.function = function(x) x)

  output$overviewTable <- renderTable({
    myTable <- data.frame(
      `Parameter` = c("\\(\\text{Statistical power}\\)",
                      "\\(\\text{Number of clusters}\\)",
                      "\\(\\text{Cluster size}\\)",
                      "\\(\\text{Family-wise false positive rate}\\)",
                      "\\(\\text{Effect for }Y_1\\)",
                      "\\(\\text{Effect for }Y_2\\)",
                      "\\(\\text{Total variance of }Y_1\\)",
                      "\\(\\text{Total variance of }Y_2\\)",
                      "\\(\\text{Endpoint-specific ICC for }Y_1\\)",
                      "\\(\\text{Endpoint-specific ICC for }Y_2\\)",
                      "\\(\\text{Inter-subject between-endpoint ICC}\\)",
                      "\\(\\text{Intra-subject between-endpoint ICC}\\)",
                      "\\(\\text{Treatment allocation ratio}\\)"),
      `Notation` = c("\\(\\pi\\)",
                                 "\\(K\\)",
                                 "\\(m\\)",
                                 "\\(\\alpha\\)",
                                 "\\(\\beta_1^*\\)",
                                 "\\(\\beta_2^*\\)",
                                 "\\(\\sigma_1^2\\)",
                                 "\\(\\sigma_2^2\\)",
                                 "\\(\\rho_0^{(1)}\\)",
                                 "\\(\\rho_0^{(2)}\\)",
                                 "\\(\\rho_1^{(1,2)}\\)",
                                 "\\(\\rho_2^{(1,2)}\\)",
                                 "\\(r\\)"),
      `Variable` = c("\\(\\text{power}\\)",
                                     "\\(\\text{K}\\)",
                                     "\\(\\text{m}\\)",
                                     "\\(\\text{alpha}\\)",
                                     "\\(\\text{beta1}\\)",
                                     "\\(\\text{beta2}\\)",
                                     "\\(\\text{varY1}\\)",
                                     "\\(\\text{varY2}\\)",
                                     "\\(\\text{rho01}\\)",
                                     "\\(\\text{rho02}\\)",
                                     "\\(\\text{rho1}\\)",
                                     "\\(\\text{rho2}\\)",
                                     "\\(\\text{r}\\)"),
      `Description` = c("\\(\\text{Probability of detecting a true effect under } H_A\\)",
                        "\\(\\text{Number of clusters in each treatment arm}^1\\)",
                        "\\(\\text{Number of individuals in each cluster}\\)",
                        "\\(\\text{Probability of one or more Type I error(s)}\\)",
                        "\\(\\text{Estimated intervention effect on the first outcome }(Y_1)\\)",
                        "\\(\\text{Estimated intervention effect on the second outcome }(Y_2)\\)",
                        "\\(\\text{Total variance of the first outcome, } Y_1\\)",
                        "\\(\\text{Total variance of the second outcome, } Y_2\\)",
                        "\\(\\text{Correlation for } Y_1 \\text{ for two different individuals in the same cluster}\\)",
                        "\\(\\text{Correlation for } Y_2 \\text{ for two different individuals in the same cluster}\\)",
                        "\\(\\text{Correlation between } Y_1 \\text{ and } Y_2 \\text{ for two different individuals in the same cluster}\\)",
                        "\\(\\text{Correlation between } Y_1 \\text{ and } Y_2 \\text{ for the same individual}\\)",
                        "\\(\\text{Treatment allocation ratio; } K_2 = rK_1\\)")
    )}, sanitize.text.function = function(x) x)



}

# Run the application
shinyApp(ui = ui, server = server)
