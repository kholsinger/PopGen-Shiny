library(plotly)

initial_p <- 0.5
initial_N <- 100
initial_n_gen <- 100
initial_n_pops <- 5

## from https://plot.ly/r/cumulative-animations/
##
accumulate_by <- function(dat, var) {
  var <- lazyeval::f_eval(var, dat)
  lvls <- plotly:::getLevels(var)
  dats <- lapply(seq_along(lvls), function(x) {
    cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
  })
  dplyr::bind_rows(dats)
}

## simulate one generation of drift
##
one_generation <- function(p, N) {
  k <- rbinom(1, N, p)
  return(k/N)
}

## return drift simulation for one population
##
drift <- function(p0, N, n_gen) {
  p <- numeric(n_gen + 1)
  p[1] <- p0
  for (i in 1:n_gen) {
    p[i+1] <- one_generation(p[i], N)
  }
  return(p)
}


df <- data.frame(N = NULL,
                 p = NULL,
                 pop = NULL)
for (i in 1:n_pops) {
  N <- seq(from = 0, to = initial_n_gen, by = 1)
  p <- drift(initial_p, initial_N, initial_n_gen)
  pop <- rep(paste("Pop", i, sep=""), length(N))
  tmp <- data.frame(N = N,
                    p = p,
                    pop = pop)
  df <- rbind(df, tmp)
}

d <- df %>%
  accumulate_by(~N)

p_plot <- d %>%
  plot_ly(
    x = ~N,
    y = ~p,
    split = ~pop,
    frame = ~frame,
    mode = "lines",
    type = "scatter",
    line = list(simplyfy = FALSE),
    showlegend = FALSE) %>%
  layout(
    yaxis = list(range = c(0.0, 1.0))) %>%
  animation_opts(
    frame = 100,
    transition = 0,
    redraw = FALSE
  ) %>%
  animation_slider(
    hide = TRUE
  )
print(p_plot)
