library(shiny)
library(plotly)

## adjustable via sliders
##
initial_p <- 0.5
initial_N <- 100
initial_n_gen <- 100
initial_n_pops <- 5

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
                  value = initial_n_gen),
      sliderInput("n_pops",
                  "Number of populations",
                  min = 1,
                  max = 10,
                  value = initial_n_pops)
      ),

    ## Show a plot of the simulated change in allele frequencies
    ##
    mainPanel(
      p("Notes explaining the principles of genetic drift are available at:",
        uiOutput("darwin")),
      p("The sliders to the left allow you to select a different initial allele frequency, population size, number of generations, and number of populations. Each line represents the history of allele frequency change in one population. All populations begin with an identical allele frequency. Each time you change one of the sliders, you'll get a new set of simulation results. If you hit \"Play\" without changing a slider, you'll get a duplicate of the plot you just saw."),
      plotlyOutput("allele_frequency_plot"),
      hr(),
      p("Source code for this and other Shiny applications is available at:",
        uiOutput("github"))
    )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  url_1 <- a("http://darwin.eeb.uconn.edu/eeb348-notes/drift.pdf",
           href="http://darwin.eeb.uconn.edu/eeb348-notes/drift.pdf")
  output$darwin <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
           href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })

  output$allele_frequency_plot <- renderPlotly({
    ## generate allele frequencies
    ##
    df <- data.frame(N = NULL,
                     p = NULL,
                     pop = NULL)
    for (i in 1:input$n_pops) {
      t <- seq(from = 0, to = input$n_gen, by = 1)
      p <- drift(input$p, input$N, input$n_gen)
      pop <- rep(paste("Pop", i, sep=""), length(t))
      tmp <- data.frame(t = t,
                        p = p,
                        pop = pop)
      df <- rbind(df, tmp)
    }
    ## construct data frame for plot
    ##
    d <- df %>%
      accumulate_by(~t)
    ## plot it
    ##
    p_plot <- d %>%
      plot_ly(
        x = ~t,
        y = ~p,
        split = ~pop,
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

