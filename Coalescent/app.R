library(ggplot2)
library(shiny)
library(cowplot)
library(plotly)
library(coala)
library(ape)

rm(list=ls())

tree_env <- new.env()
parent.env(tree_env)

## number of haplotypes to select
##
initial_N <- 10

## Define UI
##
ui <- fluidPage(
  titlePanel("The genealogy of the coalescent at one locus in one population"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("N",
                  "Number of alleles",
                  min = 1,
                  max = 100,
                  value = initial_N),
      actionButton("go", "Go")
    ),
    mainPanel(
      p("Notes explaining principles of the coalescent are available at:"),
      uiOutput("darwin"),
      p("Use the slider to select the number of alleles in the sample."),
      p("Hit \"Go\" to run the simulation."),
      p("Notice that both the topology and depth of the tree (distance from the any tip to the base) differ from simulation to simulation, even if you don't change the number of alleles sampled."),
      h4("Coalescent tree"),
      p("The coalescent tree is displayed in such a way that the distance from the tips of the tree to the base is always the same regardeless of the actual tree depth."),
      plotOutput("plot"),
      h4("Interpreting the tree depth"),
      textOutput("text"),
      hr(),
      p("Source code for this and other Shiny applications is available at:"),
      uiOutput("github")
    )
  )
)

## Define server logic
##
server <- function(input, output) {
  url_1 <- a("http://darwin.eeb.uconn.edu/eeb348/lecture-notes/coalescent.pdf",
           href="http://darwin.eeb.uconn.edu/eeb348/lecture-notes/coalescent.pdf")
  output$darwin <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
           href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })

  get_tree <- eventReactive(input$go, {
    haplos <- coal_model(sample_size = input$N,
                         loci_number = 1) +
      sumstat_trees()
    tree <- simulate(haplos)
    ape_tree <- read.tree(text = tree$trees[[1]])
    tree_env$tree_depth <- max(branching.times(ape_tree))
    title <- paste("Tree depth: ", round(tree_env$tree_depth, 3))
    if (input$N <= 20) {
      plot(ape_tree, main = title)
    } else {
      plot(ape_tree, show.tip.label = FALSE, main = title)
    }
  })

  output$plot <- renderPlot({
    get_tree()
  })

  get_text <- eventReactive(input$go, {
    paste("The actual tree depth of ",
          round(tree_env$tree_depth, 3),
          " is measured in units 4Ne generations ",
          "(assuming that the population is diploid). ",
          "Thus, for Ne = 100, the number of generations from each tip ",
          "to the base of the tree is ",
          round(4*100*tree_env$tree_depth, 3),
          " generations.\n",
          "The fractional number of generations arises ",
          "because the simulation model uses the continuous ",
          "version of the coalescent process that is described ",
          "in the notes.\n",
          "Note: The tree depth statistics reported here are rounded to ",
          "3 decimal places. If you multiply the reported tree depth by ",
          "400, you are not likely to get the number reported here",
          "(although it will be close).",
          sep = "")
  })

  output$text <- renderText({
    get_text()
  })
}

## Run the application
##
shinyApp(ui = ui, server = server)

