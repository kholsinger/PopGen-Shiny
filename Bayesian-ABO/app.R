library(shiny)
library(bayesplot)
library(R2jags)
library(ggplot2)

rm(list=ls())

initial_counts <- c(862, 131, 365, 702)
initial_total <- sum(initial_counts)

result <- new.env()
parent.env(result)

run_jags <- function(x, N) {
  jags.data <- c("x", "N")
  jags.par <- c("p.a", "p.b", "p.o")

  fit <- jags(data=jags.data,
              model.file="multinomial.jags",
              parameters.to.save=jags.par,
              n.burnin=1000,
              n.iter=2000,
              n.thin=1,
              n.chains=5)
  return(fit)
}

get_parameter <- function(fit, par) {
  x <- fit$BUGSoutput$sims.list[[par]]
  x_bar <- mean(x)
  limit <- quantile(x, c(0.025, 0.975))
  return(data.frame(x_bar = x_bar,
                    lo = limit[1],
                    hi = limit[2]))
}

estimate <- function(na, nab, nb, no) {
  x <- c(na, nab, nb, no)
  N <- sum(x)
  fit <- run_jags(x, N)
  result$fit <- fit
  p.a <- get_parameter(fit, "p.a")
  p.b <- get_parameter(fit, "p.b")
  p.o <- get_parameter(fit, "p.o")
  dat <- data.frame(Statistic = c("Mean", "2.5%", "97.5%"),
                    p.a = c(round(p.a$x_bar, 3),
                      round(p.a$lo, 3),
                      round(p.a$hi, 3)),
                    p.b = c(round(p.b$x_bar, 3),
                      round(p.b$lo, 3),
                      round(p.b$hi, 3)),
                    p.o = c(round(p.o$x_bar, 3),
                      round(p.o$lo, 3),
                      round(p.o$hi, 3)))
  return(dat)
}

## Define UI
##
ui <- fluidPage(
  titlePanel("Estimating allele frequencies in the ABO blood group system (Bayesian version)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("na",
                  "Number of type A",
                  min=1,
                  max=1000,
                  value=initial_counts[1]),
      sliderInput("nab",
                  "Number of type AB",
                  min=1,
                  max=1000,
                  value=initial_counts[2]),
      sliderInput("nb",
                  "Number of type B",
                  min=1,
                  max=1000,
                  value=initial_counts[3]),
      sliderInput("no",
                  "Number of type O",
                  min=1,
                  max=1000,
                  value=initial_counts[4]),
      actionButton("go", "Go")
    ),
    mainPanel(
      p("Notes explaining the data and the algorithm are available at:"),
      uiOutput("darwin"),
      p("Calculations will begin when you press \"Go\"."),
      h2("Phenotype counts"),
      fluidRow(
        column(4,
               dataTableOutput("phenos")
               )
        ),
      h2("Estimated allele frequencies"),
      fluidRow(
        column(4,
               dataTableOutput("alleles")
               )
        ),
      h2("Posterior distributions of allele frequencies"),
      plotOutput("plot"),
      p("Source code for this and other Shiny applications is available at:"),
      uiOutput("github")
    )
  )
)

## Define server logic
##
server <- function(input, output, session) {
  url_1 <- a("http://darwin.eeb.uconn.edu/eeb348/lecture-notes/hardy-weinberg.pdf",
           href="http://darwin.eeb.uconn.edu/eeb348/lecture-notes/hardy-weinberg.pdf")
  output$darwin <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
           href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })

  get_phenos <- eventReactive(input$go,
                              return(data.frame(A = input$na,
                                                AB = input$nab,
                                                B = input$nb,
                                                O = input$no)))
  output$phenos <- renderDataTable(get_phenos(),
                                   options=list("paging"=FALSE,
                                     "ordering"=FALSE,
                                     "info"=FALSE,
                                     "searching"=FALSE))

  get_alleles <- eventReactive(input$go,
                               estimate(input$na,
                                        input$nab,
                                        input$nb,
                                        input$no))
  output$alleles <- renderDataTable(get_alleles(),
                                    options=list("paging"=FALSE,
                                      "ordering"=FALSE,
                                      "info"=FALSE,
                                      "searching"=FALSE))

  get_plot <- eventReactive(input$go, {
    outfile <- tempfile(file = ".png")
    png(outfile, width = 800, height = 400)
    p <- mcmc_hist(result$fit$BUGSoutput$sims.array,
                   pars = c("p.a","p.b","p.o"),
                   binwidth = 0.01) +
      xlim(0.0, 1.0)
    print(p)
    dev.off()
    return(list(src = outfile,
                alt = "Posterior distributions of allele frequencies"))
  })

  output$plot <- renderImage({
    get_plot()
  })
}

## Run the application
##
shinyApp(ui = ui, server = server)

