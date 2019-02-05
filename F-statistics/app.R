library(shiny)
library(pegas)

rm(list=ls())

## adjustable via sliders
##
initial_fst <- 0.1
initial_fis <- 0.0
initial_p <- 0.5
## fixed for all examples
##
n <- 25
n_pops <- 5

column_names <- function(n_pops) {
  cnames <- character(n_pops)
  for (i in 1:n_pops) {
    cnames[i] <- paste("Pop-", i, sep="")
  }
  return(cnames)
}

sample_genos <- function(fst, fis, p) {
  a <- ((1.0 - fst)/fst)*p
  b <- ((1.0 - fst)/fst)*(1-p)
  p_pop <- rbeta(n_pops, a, b)
  x <- matrix(nrow = 3, ncol = n_pops)
  x[1,] <- p_pop^2 + fis*p*(1.0 - p)
  x[2,] <- 2.0*p_pop*(1.0 - p_pop)*(1.0 - fis)
  x[3,] <- (1.0 - p_pop)^2 + fis*p*(1.0 - p)
  y <- matrix(nrow = 3, ncol = n_pops)
  for (i in 1:n_pops) {
    tmp <- rmultinom(1, n, x[,i])
    for (j in 1:3) {
      y[j,i] <- tmp[j]
    }
  }
  dat <- as.data.frame(y)
  colnames(dat) <- column_names(n_pops)
  dat$Genotype <- c("AA", "AB", "BB")
  dat <- dat[c("Genotype", colnames(dat))]
  dat$Genotype.1 <- NULL
  return(dat)
}

get_p_bar <- function(dat) {
  num <- 2*sum(dat[1,-1]) + sum(dat[2,-1])
  den <- 2*(sum(dat[1,-1]) + sum(dat[2,-1]) + sum(dat[3,-1]))
  return(num/den)
}

diversity <- function(p, n) {
  h <- (n/(n - 1.0))*2.0*p*(1.0 - p)
  return(h)
}

get_p <- function(x) {
  p <- (2.0*x[1] + x[2])/(2.0*sum(x))
  return(p)
}

gst <- function(dat) {
  h_i <- numeric(n_pops)
  h_s <- numeric(n_pops)
  p_bar <- get_p_bar(dat)
  n_tot <- n_pops*n
  h_t <- diversity(p_bar, n_tot)
  for (i in 2:(n_pops+1)) {
    p_i <- get_p(dat[,i])
    h_s[i] <- diversity(p_i, n)
    h_i[i] <- dat[2,i]/sum(dat[,i])
  }
  f_is <- 1.0 - mean(h_i)/mean(h_s)
  f_st <- 1.0 - mean(h_s)/h_t
  f_it <- 1.0 - mean(h_i)/h_t
  return(list(f_is = f_is,
              f_st = f_st,
              f_it = f_it))
}

theta <- function(dat) {
  pop <- character(0)
  x <- numeric(0)
  for (i in 1:n_pops) {
    pop <- c(pop, rep(paste("Pop-", i, sep=""), n))
    ## dat includes a genotype indicator in column 1
    ##
    ct <- dat[,i+1]
    x <- c(x, rep("a/a", ct[1]))
    x <- c(x, rep("a/b", ct[2]))
    x <- c(x, rep("b/b", ct[3]))
  }
  alleles <- data.frame(Population = pop,
                        Locus = x)
  weir <- Fst(as.loci(alleles))
  return(list(f_is = weir[1,"Fis"],
              f_st = weir[1,"Fst"],
              f_it = weir[1,"Fit"]))
}

get_fstats <- function(dat) {
  nei <- gst(dat)
  weir <- theta(dat)
  fstats <- data.frame(Estimate =c("Nei", "Weir & Cockerham"),
                       Fis = c(round(nei$f_is, 3),
                               round(weir$f_is, 3)),
                       Fst = c(round(nei$f_st, 3),
                               round(weir$f_st, 3)),
                       Fit = c(round(nei$f_it, 3),
                               round(weir$f_it, 3)))
  return(fstats)
}

## Define UI
##
ui <- fluidPage(
  titlePanel("Genetics of geographically structured populations"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("fst",
                  "Parametric Fst:",
                  min = 0.0,
                  max = 1.0,
                  value = initial_fst),
      sliderInput("fis",
                  "Parametric Fis:",
                  min = 0.0,
                  max = 1.0,
                  value = initial_fis),
      sliderInput("p",
                  "p_bar:",
                  min = 0.0,
                  max = 1.0,
                  value = initial_p),
      actionButton("go", "Go")
    ),
    mainPanel(
      p("Notes explaining the data and the algorithm are available at:"),
      uiOutput("darwin"),
      p("Notice that each time you hit \"Go\" you'll get a different set of genotype counts, even if you don't change the Parametric Fst. That's both because there's sampling of genotypes from the underlying population genotype frequencies and because the population genotype frequencies are themselves derived from allele frequencies that are sampled from the parametric allele frequency distribution among populations."),
      h4("Genotype frequencies in populations"),
      fluidRow(
        column(n_pops + 1,
               dataTableOutput("genotypes")
               )
      ),
      h4("F-statistics"),
      fluidRow(
        column(3,
               dataTableOutput("fstats")
        )
      ),
      p("Notice that the Weir & Cockerham estimates for Fst are consistently smaller than the Nei estimates. That's because the Weir & Cockerham estimates are properly accounting for both sources of variation (sampling of populations",
        strong(em("and")), "sampling of genotypes within populations). The Nei estimates are accounting only for the within population sampling."),
      p("If you hit \"Go\" a number of times, you may also notice that the Weir & Cockerham estimates are closer (on average) to the parametric value you've chosen than the Nei estimates."),
      p("Implementation note: I hand-coded the Nei estimates. I used Fst from pegas to get the Weir & Cockerham estimates."),
      hr(),
      p("Source code for this and other Shiny applications is available at:"),
      uiOutput("github")
    )
  )
)

## Define server logic
##
server <- function(input, output) {
  url_1 <- a("http://darwin.eeb.uconn.edu/eeb348/lecture-notes/genetic-structure.pdf",
           href="http://darwin.eeb.uconn.edu/eeb348/lecture-notes/genetic-structure.pdf")
  output$darwin <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
           href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })

  current_seed <- as.numeric(Sys.time())
  ## get genotype counts in sampled populations
  ##
  set.seed(current_seed)
  get_genos <- eventReactive(input$go,
                             sample_genos(input$fst,
                                          input$fis,
                                          input$p))
  output$genotypes <- renderDataTable(get_genos(),
                                      options=list("paging"=FALSE,
                                                   "ordering"=FALSE,
                                                   "info"=FALSE,
                                                   "searching"=FALSE))
  ## get F-statistics
  ##
  set.seed(current_seed)
  output$fstats <- renderDataTable(get_fstats(get_genos()),
                                   options=list("paging"=FALSE,
                                                "ordering"=FALSE,
                                                "info"=FALSE,
                                                "searching"=FALSE))
}

## Run the application
##
shinyApp(ui = ui, server = server)

