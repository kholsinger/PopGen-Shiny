library(R2jags)
library(ggplot2)
library(bayesplot)

x <- c(862, 131, 365, 702)
N <- sum(x)

jags.data <- c("x", "N")
jags.par <- c("p.a", "p.b", "p.o")

fit <- jags(data=jags.data,
            model.file="multinomial.jags",
            parameters.to.save=jags.par,
            n.burnin=1000,
            n.iter=2000,
            n.thin=1,
            n.chains=5)

opt.old <- options(width=180)
print(fit, digits=3)
options(opt.old)

p <- mcmc_hist(fit$BUGSoutput$sims.array,
               pars = c("p.a","p.b","p.o"),
               binwidth = 0.01) +
  xlim(0.0, 1.0)
print(p)
