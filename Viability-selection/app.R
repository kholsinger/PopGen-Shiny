library(ggplot2)
library(shiny)
library(cowplot)
library(plotly)

rm(list=ls())

## adjustable via sliders
##
initial_w11 <- 0.6
initial_w12 <- 0.9
initial_w22 <- 0.45
initial_n_gen <- 25

## fixed starting points
##
p0 <- c(0.2, 0.5, 0.8)

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

## dynamics of viability selection in one population
##
viability_selection <- function(p0, w11, w12, w22, n_gen) {
  p <- numeric(n_gen + 1)
  p[1] <- p0
  for (i in 2:(n_gen + 1)) {
    q <- 1.0 - p[i-1]
    w_bar <- p[i-1]^2*w11 + 2*p[i-1]*q*w12 + q^2*w22
    p[i] <- (p[i-1]^2*w11 + p[i-1]*q*w12)/w_bar
  }
  return(p)
}

## Define UI
##
ui <- fluidPage(
  titlePanel("Viability selection at one locus with two alleles"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("w11",
                  "w11",
                  min = 0.0,
                  max = 1.0,
                  value = initial_w11),
      sliderInput("w12",
                  "w12",
                  min = 0.0,
                  max = 1.0,
                  value = initial_w12),
      sliderInput("w22",
                  "w22",
                  min = 0.0,
                  max = 1.0,
                  value = initial_w22),
      sliderInput("n_gen",
                  "Number of generations",
                  min = 1,
                  max = 100,
                  value = initial_n_gen)
    ),
    mainPanel(
      p("Notes explaining the data and calculations are available at:"),
      uiOutput("darwin"),
      p("If you select fitnesses that are either overdominant or underdominant, a dashed line from the x-axis to the mean fitness curve will show the location of the polymorphic equilibrium."),
      h4("Mean fitness"),
      plotOutput("plot"),
      h4("Allele frequency dynamics"),
      p("Hit \"Play\" to see how allele frequencies change over time from three different starting points: p = 0.2, p = 0.5, and p = 0.8."),
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
  url_1 <- a("http://darwin.eeb.uconn.edu/eeb348/lecture-notes/selection.pdf",
           href="http://darwin.eeb.uconn.edu/eeb348/lecture-notes/selection.pdf")
  output$darwin <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
           href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })

  output$plot <- renderImage({
    p <- seq(from = 0.0,
             to = 1.0,
             by = 0.001)
    w_bar <- p^2*input$w11 + 2.0*p*(1.0 - p)*input$w12 + (1.0 - p)^2*input$w22
    w11 <- input$w11
    w12 <- input$w12
    w22 <- input$w22
    if (((w12 < w11) && (w12 < w22)) || ((w12 > w11) && (w12 > w22))) {
      s1 <- 1 - w11/w12
      s2 <- 1 - w22/w12
      p_eq <- s2/(s1 + s2)
      w_bar_eq <- p_eq^2*w11 + 2.0*p_eq*(1.0 - p_eq)*w12 + (1.0 - p_eq)^2*w22
    } else {
      p_eq <- NA
    }
    outfile <- tempfile(fileext = ".png")
    png(outfile, width = 400, height = 400)
    for.plot <- data.frame(p = p,
                           w_bar = w_bar)
    w_plot <- ggplot(for.plot, aes(x = p, y = w_bar)) + geom_line() +
      ylab(expression(bar(w)))
    if (!is.na(p_eq)) {
      w_plot <- w_plot + geom_segment(aes(x = p_eq, y = 0,
                                          xend = p_eq, yend = w_bar_eq),
                                      linetype = "dashed") +
        annotate("text", x = p_eq, y = w_bar_eq + 0.04,
                 label = "hat(p)",
                 parse = TRUE, size = 5)
    }
    print(w_plot)
    dev.off()
    return(list(src = outfile,
                alt = "A graph showing w-bar as a function of p"))
  }, deleteFile = TRUE)

  output$dynamic_plot <- renderPlotly({
    ## generate allele frequencies
    ##
    df <- data.frame(t = NULL,
                     p = NULL,
                     pop = NULL)
    for (i in 1:length(p0)) {
      t <- seq(from = 0, to = input$n_gen, by = 1)
      p <- viability_selection(p0[i],
                               input$w11,
                               input$w12,
                               input$w22,
                               input$n_gen)
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

## Run the application
##
shinyApp(ui = ui, server = server)

