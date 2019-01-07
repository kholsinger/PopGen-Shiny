library(shiny)
library(reshape2)
library(ggplot2)
library(cowplot)

rm(list=ls())

init_p <- 0.5
init_AA <- 100
init_Aa <- 50
init_aa <- 0
init_ve <- 25^2
init_sample_size <- log10(1000)

variances <- new.env()
## initialize to nonsense values
##
variances$va <- -999.0
variances$vg <- -999.0
variances$ve <- -999.0
variances$x <- c(-999.0, -999.0, -999.0)
parent.env(variances)


get_variances <- function(p, x11, x12, x22, ve) {
  x_bar <- p^2*x11 + 2.0*p*(1.0 - p)*x12 + (1.0 - p)^2*x22
  a1 <- p*x11 + (1.0 - p)*x12 - x_bar/2.0
  a2 <- p*x12 + (1.0 - p)*x22 - x_bar/2.0
  vg <- p^2*(x11 - x_bar)^2 + 2.0*p*(1.0 - p)*(x12 - x_bar)^2 +
        (1.0 - p)^2*(x22 - x_bar)^2
  va <- p^2*(2*a1 - x_bar)^2 + 2.0*p*(1.0 - p)*(a1+a2 - x_bar)^2 +
        (1.0 - p)^2*(2*a2 - x_bar)^2
  vd <- p^2*(2*a1 - x11)^2 + 2.0*p*(1.0 - p)*(a1+a2 - x12)^2 +
        (1.0 - p)^2*(2*a2 - x22)^2
  variances$vg <- vg
  variances$va <- va
  variances$ve <- ve
  variances$x <- c(x11, x12, x22)
  dat <- data.frame(Component = c("Vp", "Vg", "Va", "Vd", "Ve"),
                    variance = c(round(vg + ve, 3),
                                 round(vg, 3),
                                 round(va, 3),
                                 round(vd, 3),
                                 round(ve, 3)))
  return(dat)
}

## genotypes:
##
## AA - 3
## Aa - 2
## aa - 1
##
## gamete returned
##
## A - 1
## a - 0
##
get_gamete <- function(x) {
  if (x == 3) {
    retval <- 1
  } else if (x == 1) {
    retval <- 0
  } else {
    if (runif(1) < 0.5) {
      retval <- 1
    } else {
      retval <- 0
    }
  }
  return(retval)
}

get_offspring <- function(x, y) {
  a <- get_gamete(x)
  b <- get_gamete(y)
  ## + 1 because genotypes go from 1..3
  ##
  return(a + b + 1)
}


## from Chang, Winston. R Graphics Cookbook. O'Reilly Media.
##
## Given a model, predict values of yvar from xvar
## This supports one predictor and one predicted variable
##
## xrange: If NULL, determine the x range from the model object. If a vector
##    with two numbers, use those as the min and max of the prediction range.
##
## samples: Number of samples across the x range.
##
## ...: Further arguments to be passed to predict()
##
predictvals <- function(model, xvar, yvar, xrange = NULL, samples = 100, ...) {  ## If xrange isn't passed in, determine xrange from the models.
  ## Different ways of extracting the x range, depending on model type
  if (is.null(xrange)) {
    if (any(class(model) %in% c("lm", "glm"))) {
      xrange <- range( model $ model[[ xvar]])
    } else if (any(class( model) %in% "loess")) {
      xrange <- range(model$x)
    }
  }
  newdata <- data.frame(x = seq(xrange[1], xrange[2], length.out = samples))
  names(newdata) <- xvar
  newdata[[yvar]] <- predict(model, newdata = newdata, ...)
  return(newdata)
}

## Define UI
##
ui <- fluidPage(
  titlePanel("Resemblance between parents and offspring"),
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
                  value = init_aa),
      sliderInput("ve",
                  "Environmental variance",
                  min = 1.0,
                  max = 10000.0,
                  value = init_ve),
      sliderInput("sample_size",
                  "log10(Sample size) for parent-offspring regression",
                  min = log10(10),
                  max = log10(100000),
                  value = init_sample_size)
    ),
    mainPanel(
      p("This application illustrates the resemblance between parents and offspring. It uses the allele frequency and genotypic values to calculate the additive and dominance variance. The phenotypic variance is the sum of those and the environmental variance. You'll find more information in the notes at"),
      uiOutput("darwin"),
      h2("Phenotypic distribution"),
      plotOutput("phenotypes"),
      h2("Variance components"),
      fluidRow(
        column(2,
               dataTableOutput("variances")
               )
      ),
      h2("Heritability"),
      uiOutput("narrow"),
      uiOutput("broad"),
      h2("Parent-offspring regression"),
      p("To illustrate the relationship between parents and offspring, we"),
      HTML("<ol>
            <li>Generate females and males at random by (a) picking a genotypic values at random based on the underlying genotype frequencies and (b) assigning by at random from a normal distribution with a mean given by the genotypic value and a variance give by Ve.</li>
            <li>Make male-female pairs at random.</li>
            <li>Produce 1 offspring per pair using Mendel's rules.</li>
            <li>Assign the offspring a phenotype based on its genotype and Ve.</li>
            <li>Plot the offspring on the y-axis vs. the mid-parent value on the x-axis and calculate the regression. The slope of the regression should be approximately equal to the heritability calculated above</li>
            </ol>"),
      plotOutput("parent_offspring"),
      hr(),
      p("Source code for this and other Shiny applications is available at:"),
      uiOutput("github")
    )
  )
)

## Define server logic
##
server <- function(input, output, session) {
  url_1 <- a("http://darwin.eeb.uconn.edu/eeb348-notes/quant-resemblance.pdf",
             href="http://darwin.eeb.uconn.edu/eeb348-notes/quant-resemblance.pdf")
  output$darwin <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
           href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })

  output$phenotypes <- renderImage({
    minimum <- min(input$AA, input$Aa, input$aa) - 4*sqrt(input$ve)
    maximum <- max(input$AA, input$Aa, input$aa) + 4*sqrt(input$ve)
    x <- seq(from = minimum,
             to = maximum,
             by = 0.01)
    AA <- dnorm(x, mean = input$AA, sd = sqrt(input$ve))
    Aa <- dnorm(x, mean = input$Aa, sd = sqrt(input$ve))
    aa <- dnorm(x, mean = input$aa, sd = sqrt(input$ve))
    p <- input$p
    Population <- p^2*AA + 2.0*p*(1.0 - p)*Aa + (1.0 - p)^2*aa
    wide.for.plot <- data.frame(x = x,
                                AA = p^2*AA,
                                Aa = 2.0*p*(1.0 - p)*Aa,
                                aa = (1.0 - p)^2*aa,
                                Population = Population)
    for.plot <- melt(wide.for.plot,
                     id.vars = "x",
                     variable.name = "Genotype",
                     value.name = "y")
    outfile <- tempfile(fileext = ".png")
    png(outfile, width = 800, height = 400)
    pop_plot <- ggplot(for.plot, aes(x = x, y = y, color = Genotype)) +
      geom_line() +
      scale_fill_discrete(limits = c("AA", "Aa", "aa", "Population")) +
      labs(x = "Phenotype") +
      theme(legend.title = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank()) +
      scale_y_continuous(breaks = NULL)
    print(pop_plot)
    dev.off()
    return(list(src = outfile,
                alt = "A graph showing phenotype distributions"))
  })

  output$variances <- renderDataTable(get_variances(input$p,
                                                    input$AA,
                                                    input$Aa,
                                                    input$aa,
                                                    input$ve),
                                      options=list("paging"=FALSE,
                                                   "ordering"=FALSE,
                                                   "info"=FALSE,
                                                   "searching"=FALSE))

  output$narrow <- renderText({
    variances <- get_variances(input$p,
                               input$AA,
                               input$Aa,
                               input$aa,
                               input$ve)
    va <- variances$variance[variances$Component == "Va"]
    vg <- variances$variance[variances$Component == "Vg"]
    ve <- variances$variance[variances$Component == "Ve"]
    HTML(paste("h<sup>2</sup><sub>N</sub> = "),
          round(va/(vg + ve), 3))
  })

  output$broad <- renderText({
    variances <- get_variances(input$p,
                               input$AA,
                               input$Aa,
                               input$aa,
                               input$ve)
    vg <- variances$variance[variances$Component == "Vg"]
    ve <- variances$variance[variances$Component == "Ve"]
    HTML(paste("h<sup>2</sup><sub>B</sub> = ",
               round(vg/(vg + ve), 3)))
  })

  output$parent_offspring <- renderImage({
    p <- input$p
    ve <- input$ve
    sample_size <- 10^input$sample_size
    genotypes <- c(p^2, 2.0*p*(1.0 - p), (1.0 - p)^2)
    female <- sample(3, size = sample_size, prob = genotypes, replace = TRUE)
    male <- sample(3, size = sample_size, prob = genotypes, replace = TRUE)
    mid <- numeric(sample_size)
    off <- numeric(sample_size)
    for (i in seq(from = 1, to = sample_size)) {
      f_pheno <- rnorm(1, mean = variances$x[female[i]], sd = sqrt(ve))
      m_pheno <- rnorm(1, mean = variances$x[male[i]], sd = sqrt(ve))
      mid[i] <- (f_pheno + m_pheno)/2.0
      off_geno <- get_offspring(female[i], male[i])
      off[i] <- rnorm(1, mean = variances$x[off_geno], sd = sqrt(ve))
    }
    for.plot <- data.frame(x = mid,
                           y = off,
                           x_center = mid - mean(mid),
                           y_center = off - mean(off))
    ## regression on centered mid-parent and offspring values
    ## forces regression line through mean of each
    ##
    model <- lm(y_center ~ x_center, data = for.plot)
    slope <- summary(model)$coefficients["x_center", "Estimate"]
    intercept <- summary(model)$coefficients["(Intercept)", "Estimate"]
    ## y-intercept is on shifted scale. Must add y offset (mean(mid)) and
    ## project back to original x-axis
    ##
    intercept <- mean(mid) - mean(off)*slope
    label <- paste("h^2 == ",
                   round(slope, 3))
    p <- ggplot(for.plot, aes(x = x, y = y)) +
      geom_point() +
        geom_abline(slope = slope, intercept = intercept, color="blue") +
          geom_vline(xintercept = mean(mid), linetype = "dashed") +
            geom_hline(yintercept = mean(off), linetype = "dashed") +
              annotate("text", label = label,
                       x = min(for.plot$x) + 5, y = max(for.plot$y) - 20,
                       parse = TRUE)
    outfile <- tempfile(fileext = ".png")
    png(outfile, width = 800, height = 400)
    print(p)
    dev.off()
    return(list(src = outfile,
                alt = "A graph showing a parent-offspring regression"))
  })
}

## Run the application
##
shinyApp(ui = ui, server = server)

