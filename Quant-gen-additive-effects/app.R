library(shiny)

rm(list=ls())

init_p <- 0.5
init_AA <- 100
init_Aa <- 50
init_aa <- 0

a1_text <- function(p, x11, x12, x22) {
  x_bar <- p^2*x11 + 2.0*p*(1.0 - p)*x12 + (1.0 - p)^2*x22
  text <- paste("= ", p, "*", x11, " + ", (1.0 - p), "*", x12,
                " - ", round(x_bar/2, 3),
                " = ", round(p*x11 + (1.0 - p)*x12 - x_bar/2.0, 3),
                sep = "")
  return(text)
}

a2_text <- function(p, x11, x12, x22) {
  x_bar <- p^2*x11 + 2.0*p*(1.0 - p)*x12 + (1.0 - p)^2*x22
  text <- paste("= ", p, "*", x12, " + ", (1.0 - p), "*", x22,
                " - ", round(x_bar, 3),
                " = ", round(p*x12 + (1.0 - p)*x22 - x_bar/2.0, 3),
                sep = "")
  return(text)
}

get_genotypes <- function(p, x11, x12, x22) {
  x_bar <- p^2*x11 + 2.0*p*(1.0 - p)*x12 + (1.0 - p)^2*x22
  a1 <- p*x11 + (1.0 - p)*x12 - x_bar/2.0
  a2 <- p*x12 + (1.0 - p)*x22 - x_bar/2.0
  dat <- data.frame(Type = c("Genotypic", "Additive"),
                    AA = c(x11, round(2*a1, 3)),
                    Aa = c(x12, round(a1 + a2, 3)),
                    aa = c(x22, round(2*a2, 3)),
                    Mean = c(x_bar,
                             round(p^2*(2*a1) +
                                     2.0*p*(1.0 - p)*(a1 + a2) +
                                     (1.0 - p)^2*(2*a2), 3)))
  return(dat)
}

get_variances <- function(p, x11, x12, x22) {
  x_bar <- p^2*x11 + 2.0*p*(1.0 - p)*x12 + (1.0 - p)^2*x22
  a1 <- p*x11 + (1.0 - p)*x12 - x_bar/2.0
  a2 <- p*x12 + (1.0 - p)*x22 - x_bar/2.0
  vg <- p^2*(x11 - x_bar)^2 + 2.0*p*(1.0 - p)*(x12 - x_bar)^2 +
        (1.0 - p)^2*(x22 - x_bar)^2
  va <- p^2*(2*a1 - x_bar)^2 + 2.0*p*(1.0 - p)*(a1+a2 - x_bar)^2 +
        (1.0 - p)^2*(2*a2 - x_bar)^2
  vd <- p^2*(2*a1 - x11)^2 + 2.0*p*(1.0 - p)*(a1+a2 - x12)^2 +
        (1.0 - p)^2*(2*a2 - x22)^2
  dat <- data.frame(Component = c("Vg", "Va", "Vd"),
                    variance = c(round(vg, 3),
                                 round(va, 3),
                                 round(vd, 3)))
  return(dat)
}

## Define UI
##
ui <- fluidPage(
  titlePanel("The additive effect of alleles and partitioning genetic variance"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("p",
                  "Allele frequency",
                  min = 0.0,
                  max = 1.0,
                  value = init_p),
      sliderInput("AA",
                  "Genotypic value of AA",
                  min = 0.0,
                  max = 100.0,
                  value = init_AA),
      sliderInput("Aa",
                  "Genotypic value of Aa",
                  min = 0.0,
                  max = 100.0,
                  value = init_Aa),
      sliderInput("aa",
                  "Genotypic value of aa",
                  min = 0.0,
                  max = 100.0,
                  value = init_aa)
    ),
    mainPanel(
      p("This application illustrates calculation of the additive effect of an allele and partitioning of the genotypic variance into additive and dominance components given known genotypic values for genotypes at one locus with two alleles and an allele frequency. Notes that explain the principles behind the calculation are available at"),
      uiOutput("darwin"),
      h2("Additive effects"),
      HTML("&alpha;<sub>1</sub> = p*x<sub>11</sub> + (1-p)*x<sub>12</sub> - x_bar/2"),
      uiOutput("a1"),
      HTML("&alpha;<sub>2</sub> = p*x<sub>12</sub> + (1-p)*x<sub>22</sub> - x_bar/2"),
      uiOutput("a2"),
      h2("Genotypes"),
      fluidRow(
        column(5,
               dataTableOutput("genotypes")
               )
      ),
      h2("Variance"),
      fluidRow(
        column(2,
               dataTableOutput("variances")
               )
      ),
      hr(),
      p("Source code for this and other Shiny applications is available at:"),
      uiOutput("github")
    )
  )
)

## Define server logic
##
server <- function(input, output, session) {
  url_1 <- a("http://darwin.eeb.uconn.edu/eeb348-notes/quant-intro.pdf",
             href="http://darwin.eeb.uconn.edu/eeb348-notes/quant-intro.pdf")
  output$darwin <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
           href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })

  output$a1 <- renderUI({
    a1_text(input$p, input$AA, input$Aa, input$aa)
  })
  output$a2 <- renderUI({
    a2_text(input$p, input$AA, input$Aa, input$aa)
  })

  output$genotypes <- renderDataTable(get_genotypes(input$p,
                                                    input$AA,
                                                    input$Aa,
                                                    input$aa),
                                      options=list("paging"=FALSE,
                                                   "ordering"=FALSE,
                                                   "info"=FALSE,
                                                   "searching"=FALSE))

  output$variances <- renderDataTable(get_variances(input$p,
                                                    input$AA,
                                                    input$Aa,
                                                    input$aa),
                                      options=list("paging"=FALSE,
                                                   "ordering"=FALSE,
                                                   "info"=FALSE,
                                                   "searching"=FALSE))

  observeEvent(input$exit, stopApp())
}

## Run the application
##
shinyApp(ui = ui, server = server)

