library(ggplot2)
library(shiny)
library(cowplot)
library(plotly)

rm(list=ls())

## adjustable via sliders
##
initial_p <- 0.5
initial_sigma <- 0.5
n_gen <- 10

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

## heterozygote frequencies with partial selfing
##
inbreeding <- function(p, sigma, n_gen) {
  h <- numeric(n_gen + 1)
  x <- numeric(3)
  x[1] <- p^2
  x[2] <- 2*p*(1.0 - p)
  x[3] <- (1.0 - p)^2
  h[1] <- x[2]
  for (i in 2:(n_gen + 1)) {
    y <- numeric(3)
    y[1] <- p^2*(1.0 - sigma) + (x[1] + x[2]/4.0)*sigma
    y[2] <- 2*p*(1.0 - p)*(1.0 - sigma) + (x[2]/2.0)*sigma
    y[3] <- (1.0 - p)^2*(1.0 - sigma) + (x[3] + x[2]/4.0)*sigma
    x <- y
    h[i] <- x[2]
  }
  return(h)
}

## Define UI
##
ui <- fluidPage(
  titlePanel("Decay of heterozygosity with partial selfing"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("p",
                  "p",
                  min = 0.0,
                  max = 1.0,
                  value = initial_p),
      sliderInput("sigma",
                  "sigma",
                  min = 0.0,
                  max = 1.0,
                  value = initial_sigma)
    ),
    mainPanel(
      p("Notes explaining the data and calculations are available at:"),
      uiOutput("darwin"),
      h4("Decay of heterozygosity"),
      p("p is the population allele frequency. sigma is the rate of self-fertilization."),
      p("Hit \"Play\" to see how heterozygosity changes over time."),
      plotlyOutput("dynamic_plot"),
      hr(),
      p("Source code for this and other Shiny applications is available at:"),
      uiOutput("github")
    )
  )
)

## Define server logic
##
server <- function(input, output) {
  url_1 <- a("http://darwin.eeb.uconn.edu/eeb348/lecture-notes/inbreeding.pdf",
           href="http://darwin.eeb.uconn.edu/eeb348/lecture-notes/inbreeding.pdf")
  output$darwin <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
           href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })

  output$dynamic_plot <- renderPlotly({
    ## generate heterozygosity
    ##
    df <- data.frame(t = NULL,
                     h = NULL)
    t <- seq(from = 0, to = n_gen, by = 1)
    h <- inbreeding(input$p,
                    input$sigma,
                    n_gen)
    df <- data.frame(t = t,
                     h = h)
    ## construct data frame for plot
    ##
    d <- df %>%
      accumulate_by(~t)
    ## plot it
    ##
    p_plot <- d %>%
      plot_ly(
        x = ~t,
        y = ~h,
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

## Run the application
##
shinyApp(ui = ui, server = server)

