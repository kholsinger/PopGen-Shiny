library(shiny)
library(plotly)

## adjustable via sliders
##
initial_p <- 0.5
initial_N <- 100
initial_n_gen <- 100

## from https://plot.ly/r/cumulative-animations/
##
accumulate_by <- function(dat, var) {
  var <- lazyeval::f_eval(var, dat)
  lvls <- plotly:::getLevels(var)
  dats <- lapply(seq_along(lvls), function(x) {
    cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
  })
  dplyr::bind_rows(dats)
}

## simulate one generation of drift
##
one_generation <- function(p, N) {
  k <- rbinom(1, N, p)
  return(k/N)
}

## return drift simulation for one population
##
drift <- function(p0, N, n_gen) {
  p <- numeric(n_gen + 1)
  p[1] <- p0
  for (i in 1:n_gen) {
    p[i+1] <- one_generation(p[i], N)
  }
  return(p)
}


# Define UI for application that draws a histogram
ui <- fluidPage(

  ## Application title
  ##
  titlePanel("Allele frequency changes with genetic drift"),

  ## selection initial allele frequency, population size, and number of
  ## generations for simulation
  ##
  sidebarLayout(
    sidebarPanel(
      sliderInput("p",
                  "Initial frequency of A",
                  min = 0.0,
                  max = 1.0,
                  value = initial_p),
      sliderInput("N",
                  "Population size",
                  min = 2,
                  max = 1000,
                  value = initial_N),
      sliderInput("n_gen",
                  "Number of generations",
                  min = 1,
                  max = 1000,
                  value = initial_n_gen)
      ),

    ## Show a plot of the simulated change in allele frequencies
    ##
    mainPanel(
      plotlyOutput("allele_frequency_plot")
    )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$allele_frequency_plot <- renderPlotly({
    ## generate allele frequencies
    ##
    p <- drift(input$p, input$N, input$n_gen)
    N <- seq(from = 0, to = input$n_gen, by = 1)
    ## construct data frame for plot
    ##
    df <- data.frame(p = p,
                     N = N)
    d <- df %>%
      accumulate_by(~N)
    ## plot it
    ##

    p_plot <- d %>%
      plot_ly(
        x = ~N,
        y = ~p,
        frame = ~frame,
        mode = "lines",
        type = "scatter",
        line = list(simplyfy = FALSE),
        showlegend = FALSE) %>%
      layout(
        yaxis = list(range = c(0.0, 1.0))) %>%
      animation_opts(
        frame = 100,
        transition = 0,
        redraw = FALSE
      ) %>%
      animation_slider(
        hide = TRUE
      )
  })
}

# Run the application
shinyApp(ui = ui, server = server)

