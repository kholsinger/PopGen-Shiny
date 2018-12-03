library(ggplot2)
library(shiny)
library(cowplot)

rm(list=ls())

## adjustable via sliders
##
initial_w11 <- 0.6
initial_w12 <- 0.9
initial_w22 <- 0.45

## genotype numbers (from Dobzhansky example in notes)
##
AA <- 41
AB <- 82
BB <- 27

## after selection
##
genos <- new.env()
AA_new <- AA*initial_w11
AB_new <- AB*initial_w12
BB_new <- BB*initial_w22
parent.env(genos)

## convert genotype numbers to a data frame for display
##
get_genotypes <- function(w11, w12, w22) {
  genos$AA_new <- (AA*w11)
  genos$AB_new <- (AB*w12)
  genos$BB_new <- (BB*w22)
  df <- data.frame(Time = c("Before selection", "After selection"),
                   AA = c(AA, genos$AA_new),
                   AB = c(AB, genos$AB_new),
                   BB = c(BB, genos$BB_new))
  return(df)
}

calc_p <- function(aa, ab, bb, text) {
  p <- (2*aa + ab)/(2*aa + 2*ab + 2*bb)
  tmp <- paste("Frequence of A ", text, " selection: ", round(p, 3), sep="")
  return(tmp)
}

get_w_bar <- function(w11, w12, w22) {
  print(w11)
  print(w12)
  print(w22)
  p <- seq(from = 0.0,
           to = 1.0,
           by = 0.001)
  w_bar <- p^2*w11 + 2*p*(1.0 - p)*w12 * (1.0 - p)^2*w22
  df <- data.frame(p = p,
                   w_bar = w_bar)
  return(df)
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
      actionButton("go", "Go")
    ),
    mainPanel(
      p("Notes explaining the data and calculations are available at:"),
      uiOutput("darwin"),
      p("The genotype numbers correspond to those in the Dobzhansky experiment described in the notes. You can't change them. You can change the fitnesses by moving the sliders."),
      p("This illustration currently shows the results of viability selection", em("within"), "one generation, from zygote (before selection) to adult (after selection)."),
      h4("Mean fitness"),
      plotOutput("plot"),
      p(strong("Note"), ": The genotype and allele frequency below won't update until you hit \"Go\"."),
      h4("Genotypes"),
      fluidRow(
        column(4,
               dataTableOutput("genos")
               )
      ),
      h4("Allele frequency"),
      uiOutput("before"),
      uiOutput("after"),
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

  get_genos <- eventReactive(input$go,
                             get_genotypes(input$w11,
                                           input$w12,
                                           input$w22))
  output$genos <- renderDataTable(get_genos(),
                                  options=list("paging"=FALSE,
                                               "ordering"=FALSE,
                                               "info"=FALSE,
                                               "searching"=FALSE))
  get_p_before <- eventReactive(input$go,
                                calc_p(AA,
                                       AB,
                                       BB,
                                       "before"))
  output$before <- renderText(get_p_before())
  get_p_after <- eventReactive(input$go,
                               calc_p(genos$AA_new,
                                      genos$AB_new,
                                      genos$BB_new,
                                      "after"))
  output$after<- renderText(get_p_after())
  output$plot <- renderImage({
    p <- seq(from = 0.0,
             to = 1.0,
             by = 0.001)
    w_bar <- p^2*input$w11 + 2.0*p*(1.0 - p)*input$w12 + (1.0 - p)^2*input$w22
    outfile <- tempfile(fileext = ".png")
    png(outfile, width = 400, height = 400)
    for.plot <- data.frame(p = p,
                           w_bar = w_bar)
    w_plot <- ggplot(for.plot, aes(x = p, y = w_bar)) + geom_line()
    print(w_plot)
    dev.off()
    return(list(src = outfile,
                alt = "A graph showing w-bar as a function of p"))
  }, deleteFile = TRUE)
}

## Run the application
##
shinyApp(ui = ui, server = server)

