library(shiny)

rm(list=ls())

counts <- c(862, 131, 365, 702)
total <- sum(counts)
p0 <- c(0.33, 0.33, 0.34)
iteration <- 0
n_iter <- 6

pheno_freqs <- function() {
  x <- counts/total
  return(x)
}

geno_counts <- function(p,x) {
  y <- numeric(6)
  y[1] <- counts[1]*(p[1]^2)/(p[1]^2 + 2*p[1]*p[3])
  y[2] <- counts[1]*(2*p[1]*p[3])/(p[1]^2 + 2*p[1]*p[3])
  y[3] <- counts[2]
  y[4] <- counts[3]*(p[2]^2)/(p[2]^2 + 2*p[2]*p[3])
  y[5] <- counts[3]*(2*p[2]*p[3])/(p[2]^2 + 2*p[2]*p[3])
  y[6] <- counts[4]
  stopifnot(abs(sum(y) - total) < 1)
  return(y)
}

allele_freqs <- function(y) {
  p <- numeric(3)
  p[1] <- (2*y[1] + y[2] + y[3])/(2*total)
  p[2] <- (2*y[4] + y[3] + y[5])/(2*total)
  p[3] <- (2*y[6] + y[2] + y[5])/(2*total)
  return(p)
}

html_iteration <- function(i, iterations) {
  p <- iterations[i,]
  text <- sprintf("%2d: %5.3f %5.3f %5.3f", i, p[1], p[2], p[3])
  text <- HTML(text)
  return(text)
}

run_iterations <- function(n_iter) {
  x <- pheno_freqs()
  p <- matrix(nrow=n_iter, ncol=3)
  y <- matrix(nrow=n_iter, ncol=6)
  p[1,] <- p0
  y[1,] <- geno_counts(p[1,], x)
  if (n_iter > 1) {
    for (i in 2:n_iter) {
      p[i,] <- allele_freqs(y[i-1,])
      y[i,] <- geno_counts(p[i,], x)
    }
  }
  p_table <- data.frame(Iteration=1:n_iter,
                        A=round(p[,1], 3),
                        B=round(p[,2], 3),
                        O=round(p[,3], 3))
  return(p_table)
}

## Define UI
##
ui <- fluidPage(
  titlePanel("Estimating allele frequencies with the EM algorithm"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n",
                  "Number of iterations",
                  min=1,
                  max=20,
                  value=1),
      actionButton("go", "Go")
    ),
    mainPanel(
      p("Notes explaining the data and the algorithm are available at:"),
      uiOutput("darwin"),
      h2("Phenotype counts"),
      fluidRow(
        column(4,
               dataTableOutput("phenos")
               )
      ),
      h2("Allele frequencies"),
      fluidRow(
        column(4,
               dataTableOutput("table")
               )
      )
    )
  )
)

## Define server logic
##
server <- function(input, output, session) {
  url <- a("http://darwin.eeb.uconn.edu/eeb348/lecture-notes/hardy-weinberg.pdf",
           href="http://darwin.eeb.uconn.edu/eeb348/lecture-notes/hardy-weinberg.pdf")
  output$darwin <- renderUI({
    tagList("", url)
  })
  phenos <- data.frame(A=counts[1],
                       AB=counts[2],
                       B=counts[3],
                       O=counts[4])
  output$phenos <- renderDataTable(phenos,
                                   options=list("paging"=FALSE,
                                                "ordering"=FALSE,
                                                "info"=FALSE,
                                                "searching"=FALSE))
  do_iteration <- eventReactive(input$go,
                                run_iterations(input$n))
  output$table <- renderDataTable(do_iteration(),
                                  options=list("paging"=FALSE,
                                               "ordering"=FALSE,
                                               "info"=FALSE,
                                               "searching"=FALSE))
}

## Run the application
##
shinyApp(ui = ui, server = server)

