library(ggplot2)
library(shiny)
library(cowplot)
library(coala)
## Note: May need to explicitly set Bioconductor repository as follows:
## 
## setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.8/bioc"))
## 
## where the version number (3.8 above) gets adjusted to reflect the local version of Bioconductor
## being used

library(ggtree)
library(ape)

rm(list=ls())

## sample size in each population
##
initial_N <- 5
initial_m <- log10(2.0)
initial_mu <- log10(1.0)
loci_number <- 1
locus_length <- 1000

get_char <- function(snps, haplos) {
  ret_val <- "NOT FOUND"
  for (i in 1:nrow(haplos)) {
    match <- TRUE
    for (j in 1:length(snps)) {
      ## j + 1 in haplos since letter is in the first column
      ##
      match <- match && (snps[j] == as.numeric(haplos[i, j + 1]))
    }
    if (match) {
      ret_val <- haplos[i, 1]
    }
  }
  return(ret_val)
}

simulate_tree <- function(pop_config,
                          loci_number,
                          locus_length,
                          migration_rate,
                          mutation_rate)
{
  model <- coal_model(sample_size = pop_config,
                      loci_number = loci_number,
                      loci_length = locus_length,
                      ploidy = 2) +
    feat_migration(rate = migration_rate,
                   symmetric = TRUE) +
    feat_mutation(rate = mutation_rate) +
    sumstat_seg_sites() +
    sumstat_trees()
  results <- simulate(model)

  ## translate the tree into APE format
  ##
  tree <- read.tree(text = results$trees[[1]])
  tree_depth <- max(branching.times(tree))

  snps <- results$seg_sites[[1]]$snps
  haplos <- unique(snps)
  if (nrow(haplos) < 1) {
    p <- ggplot(data.frame(x = 1,
                           y = 1,
                           label = "No polymorphism"),
                aes(x = x,
                    y = y)) +
      geom_text(aes(y = y, label = label),
                fontface = "bold",
                family = "sans",
                size = 15) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank())
    return(p)
  }
  ## add first column as character to distinguish haplotypes
  ## Note: will fail if there are more than 26 haplotypes
  ##
  haplos <- cbind(letters[seq(from = 1, to = nrow(haplos), by = 1)], haplos)
  ## associate haplotypes with tips
  ##
  haplos_char <- character(nrow(snps))
  for (i in 1:nrow(snps)) {
    haplos_char[i] <- get_char(snps[i,], haplos)
  }
  ## associate haplotypes with populations
  ##
  pops <- character(0)
  for (i in 1:length(pop_config)) {
    pops <- c(pops, rep(paste0("Pop", LETTERS[i]), pop_config[i]))
  }
  ## construct base tree
  ##
  p <- ggtree(tree) +
    theme_tree2()
  ## before adding allele and population information to it
  ##
  for_plot <- data.frame(taxa = tree$tip.label,
                         Allele = haplos_char,
                         Population = pops)
  p <- p %<+% for_plot +
    geom_tippoint(aes(color = Allele, shape = Population),
                  size = min(2.0*(50.0/sum(pop_config)), 4.0)) +
    ggtitle(paste("Tree depth: ", round(tree_depth, 3), sep="")) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold",
                                    family = "sans",
                                    size = 15))
  return(p)
}

## Define UI
##
ui <- fluidPage(
  titlePanel("The genealogy of the coalescent at one locus in two populations with migration and mutation"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("N",
                  "Sample size of individuals in each population",
                  min = 1,
                  max = 10,
                  value = initial_N),
      sliderInput("m",
                  "log10(4Nm)",
                  min = -1.0,
                  max = 2.0,
                  value = initial_m),
      sliderInput("mu",
                  "log10(4Nmu)",
                  min = -1.0,
                  max = 2.0,
                  value = initial_mu),
      actionButton("go", "Go")
    ),
    mainPanel(
      p("Notes explaining principles of the coalescent are available at:"),
      uiOutput("darwin"),
      p("Notes intended to help you interpret the results presented here are not yet available, but the simulation should be pretty self explanatory."),
      p("There are two populations designated with different symbols. An infinite sites model of mutation produces different alleles, designated with different colors. The simulation assumes that populations are diploid and that the sample size is the number of diploid individuals (so the number of alleles sampled is twice the sample size)."),
      p("Use the sliders to select the sample size within each population (assumed equal to keep things simple), the migration rate (scaled by 4 times the effective population size), and the mutation rate (also scaled by 4 times the effective population size)."),
      p("Hit \"Go\" to run the simulation."),
      p("Notice that both the topology and depth of the tree (distance from the any tip to the base) differ from simulation to simulation, even if you don't change the number of alleles sampled."),
      h4("Coalescent tree"),
      p("The coalescent tree is displayed in such a way that the distance from the tips of the tree to the base is always the same regardeless of the actual tree depth. When the mutation rate is small, you may sometimes get a \"No polymorphism\" message instead of a coalescent plot. All that means is that given the low mutation rate, this particular simulation didn't happen to produce a polymorphism. Hit \"Go\" again, and you'll get a coalescent plot, if not on the next attempt, then eventually."),
      plotOutput("plot"),
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
    pop_config <- c(input$N, input$N)
    p <- simulate_tree(pop_config,
                       loci_number,
                       locus_length,
                       10^input$m,
                       10^input$mu)
    print(p)
  })

  output$plot <- renderPlot({
    get_tree()
  })
}

## Run the application
##
shinyApp(ui = ui, server = server)
