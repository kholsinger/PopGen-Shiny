library(ggplot2)

p <- 0.5
x <- c(100, 50, 0)
sample_size <- 10000
ve <- 625

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

genotypes <- c(p^2, 2.0*p*(1.0 - p), (1.0 - p)^2)
female <- sample(3, size = sample_size, prob = genotypes, replace = TRUE)
male <- sample(3, size = sample_size, prob = genotypes, replace = TRUE)
mid <- numeric(sample_size)
off <- numeric(sample_size)
for (i in seq(from = 1, to = sample_size)) {
  f_pheno <- rnorm(1, mean = x[female[i]], sd = sqrt(ve))
  m_pheno <- rnorm(1, mean = x[male[i]], sd = sqrt(ve))
  mid[i] <- (f_pheno + m_pheno)/2.0
  off_geno <- get_offspring(female[i], male[i])
  off[i] <- rnorm(1, mean = x[off_geno], sd = sqrt(ve))
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
print(p)
