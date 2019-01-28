p <- 0.5
x <- c(0, 80, 100)
Ve <- 25^2

pop_size <- 1000

slope <- 10
intercept <- 1

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
  dat <- data.frame(Vp = vg + ve,
                    Vg = vg,
                    Va = va,
                    Vd = vd,
                    Ve = ve)
  return(dat)
}

construct_population <- function(p, x, Ve, pop_size) {
  genotypes <- c(p^2, 2.0*p*(1.0 - p), (1.0 - p)^2)
  geno_idx <- sample(3, size = pop_size, prob = genotypes, replace = TRUE)
  phenotype <- numeric(pop_size)
  for (i in 1:pop_size) {
    phenotype[i] <- rnorm(1, mean = x[geno_idx[i]], sd = sqrt(Ve))
  }
  return(data.frame(genotype = geno_idx,
                    phenotype = phenotype))
}

perform_selection <- function(pop, slope, intercept) {
  ## get heritability in full population first
  ##
  mid <- numeric(nrow(pop))
  off <- numeric(nrow(pop))
  for (i in 1:nrow(pop/2)) {
    mom <- sample(nrow(pop), size = 1)
    dad <- sample(nrow(pop), size = 1)
    mid[i] <- (pop$phenotype[mom] + pop$phenotype[dad])/2.0
    off_geno <- get_offspring(pop$genotype[mom], pop$genotype[dad])
    off[i] <- rnorm(1, mean = x[off_geno], sd = sqrt(Ve))
  }
  model <- lm(off ~ mid)
  h_2 <- summary(model)$coefficients["mid", "Estimate"]
  raw_fitness <- slope*(pop$phenotype - min(pop$phenotype)) + intercept
  ## make sure raw_fitness is non-negative
  ##
  if (min(raw_fitness) < 0) {
    raw_fitness <- raw_fitness + abs(min(raw_fitness))
  }
  p_fitness <- raw_fitness/sum(raw_fitness)
  parents <- sample(nrow(pop), prob = p_fitness, replace = TRUE)
  mid <- numeric(nrow(pop))
  off <- numeric(nrow(pop))
  for (i in 1:nrow(pop/2)) {
    mom <- sample(parents, size = 1)
    dad <- sample(parents, size = 1)
    mid[i] <- (pop$phenotype[mom] + pop$phenotype[dad])/2.0
    off_geno <- get_offspring(pop$genotype[mom], pop$genotype[dad])
    off[i] <- rnorm(1, mean = x[off_geno], sd = sqrt(Ve))
  }
  return(data.frame(before = pop$phenotype,
                    parents = parents,
                    h_2 = h_2,
                    after = mid,
                    offspring = off,
                    raw = raw_fitness,
                    geno = pop$genotype))
}

before <- construct_population(p, x, Ve, pop_size)
after <- perform_selection(before, slope, intercept)

S <- mean(after$after) - mean(after$before)
R <- mean(after$offspring) - mean(after$before)

variances <- get_variances(p, x[3], x[2], x[1], Ve)
h_2 <- variances$Va/variances$Vp

cat("Mean before selection:  ", round(mean(after$before), 3), "\n",
    "Mean after selection:   ", round(mean(after$after), 3), "\n",
    "Mean in offspring:      ", round(mean(after$offspring), 3), "\n",
    "Selection differential: ", round(S, 3), "\n",
    "Heritabiity...\n",
    "  observed:             ", round(unique(after$h_2), 3), "\n",
    "  expected:             ", round(h_2, 3), "\n",
    "Response to selection...\n",
    "  observed:             ", round(R, 3), "\n",
    "  expected (obs. h^2):  ", round(S*unique(after$h_2), 3), "\n",
    "  expected (exp. h^2):  ", round(S*h_2, 3), "\n",
    sep = ""
    )
