## The intuition behind these residuals is that if we simulate a bunch of values
## from the model and look at where each observed values falls with respect to
## the simulated values, the rank order of where the observed value falls should
## be random assuming the model is a good generative process of what we observed.
## Because where in the rank order is random, if we plotted out those rank
## orders, they'd be uniform between 1 and the number of simulations. But,
## instead of rank order we calculate the cumulative proportion (i.e., divide
## those rank orders by the number of simulations). Therefore, we should have a
## uniform(0, 1) distribution.

## Probability integral transform (PIT) / quantile residuals example for
## continuous observations

simulator <- function() {
  rnorm(n = 1000, mean = 1.4, sd = 0.7)
}

set.seed(1)

## pretend these are our observations:
observations <- simulator()

## simulate from our 'model'; pretend we've fitted it to the data
nsim <- 300 ## how many times to simulate?
sim <- matrix(nrow = length(observations), ncol = nsim)
for (j in seq(ncol(sim))) {
  sim[, j] <- simulator()
}

## get quantiles
quants <- numeric(length = nrow(sim))
for (i in seq(nrow(sim))) {
  quants[i] <- sum(observations[i] < sim[i, ]) / nsim
}

## if our model is consistent with the data, these should be uniformly
## distributed between 0 and 1
hist(quants)

## -----------------------------------------------------------------------

## let's try that with a mismatched model

make_observations <- function() {
  rnorm(n = 1000, mean = 1.4, sd = 2) # SD is mismatched
}
simulator <- function() {
  rnorm(n = 1000, mean = 1.4, sd = 0.7) # SD is mismatched
}

set.seed(1)

## pretend these are our observations:
observations <- make_observations()

## simulate from our model
nsim <- 300
sim <- matrix(nrow = length(observations), ncol = nsim)
for (j in seq(ncol(sim))) {
  sim[, j] <- simulator()
}

## get quantiles
quants <- numeric(length = nrow(sim))
for (i in seq(nrow(sim))) {
  quants[i] <- sum(observations[i] < sim[i, ]) / nsim
}

## if our model is consistent with the data, these should be uniformly
## distributed between 0 and 1
hist(quants)

## Notice how more often than expected we have observed values that fall to low
## or high in the rank order of simulated values. This happens here because the
## observations came from a process with additional dispersion (here a higher
## SD) than our simulations from our (fake) fitted model.

## -----------------------------------------------------------------------

## PIT / *randomized* quantile residuals example for integer observations:

simulator_poisson <- function() {
  rpois(n = 1000, lambda = 1.8)
}

set.seed(1)

## pretend these are our observations:
observations <- simulator_poisson()

## simulate from our model
nsim <- 300 ## how many times to simulate?
sim <- matrix(nrow = length(observations), ncol = nsim)
for (j in seq(ncol(sim))) {
  sim[, j] <- simulator_poisson()
}

## get quantiles
quants <- numeric(length = length(observations))
for (i in seq(nrow(sim))) {
  quants[i] <- sum(observations[i] < sim[i, ]) / nsim
}

## look at how these are clumped
hist(quants)
## that's because we can only simulate integer values...
## we need to add randomization for this to look uniform

## here, we'll add a line, which calculates the quantile for '<=' on top of the
## previous '<'
quants_upper <- numeric(length = length(observations))
quants_lower <- numeric(length = length(observations))
for (i in seq(nrow(sim))) {
  quants_lower[i] <- sum(observations[i] < sim[i, ]) / nsim
  quants_upper[i] <- sum(observations[i] <= sim[i, ]) / nsim
}

## now we can take a uniform draw between these bounds for every observation
## we'll do that in a vectorized fashion, but we could also loop over the
## observations again
randomized_quant <- runif(length(observations), min = quants_lower, max = quants_upper)

## now we have our uniform(0, 1) distribution of randomized quantile residuals:
hist(randomized_quant)
