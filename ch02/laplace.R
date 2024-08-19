# Modified from
# Thorson, J., and Kristensen, K. 2024. Spatio-Temporal Models for Ecologists.
# In 1st edition. Chapman and Hall/CRC, Boca Raton, FL.

library(lme4)

# prepare data (ignore everything in here):
vismba <- readRDS("ch02/vismba.rds")
samples <- data.frame("x" = vismba$gx, "y" = vismba$gy, "agb" = vismba$agb)
samples <- sf::st_as_sf(samples, coords = c("x", "y"))
grid <- sf::st_make_grid(sf::st_bbox(c(xmin = 0, xmax = 1000, ymin = 0, ymax = 500)), cellsize = c(250, 125))
grid_i <- sf::st_intersects(samples, grid)
Count_i <- tapply(samples$agb, INDEX = factor(unlist(grid_i), levels = 1:length(grid)), FUN = length)
Data <- data.frame(sf::st_coordinates(sf::st_centroid(grid)), "Count" = ifelse(is.na(Count_i), 0, Count_i))
grid_sf <- sf::st_sf(grid, Count = Count_i)
plot(grid_sf, axes = TRUE, reset = FALSE, pal = sf::sf.colors(n = 10, alpha = 0.2), breaks = seq(0, max(Data$Count), length = 11))
plot(samples, add = TRUE, pch = 20)

# fit with lme4 for comparison later:
Data$Site <- factor(1:nrow(Data))
fit_glmer <- glmer(Count ~ 1 + (1 | Site), data = Data, family = poisson(link = "log"))
summary(fit_glmer)

# function to return the joint (fixed + random) log likelihood:
joint_nll <- function(prandom, pfixed, plist, Data, random) {
  # Read in values
  phat_fixed <- relist(pfixed, plist[setdiff(names(plist), random)])
  phat_random <- relist(prandom, plist[random])
  phat <- c(phat_fixed, phat_random)
  # Data likelihood
  jnll <- 0
  for (i in 1:nrow(Data)) {
    jnll <- jnll - dpois(Data$Count[i], lambda = exp(phat$eps[i]), log = TRUE)
  }
  # Random effect distribution
  for (i in 1:nrow(Data)) {
    jnll <- jnll - dnorm(phat$eps[i], mean = phat$logmu, sd = exp(phat$logsd), log = TRUE)
  }
  jnll
}

# outer optimization function; returns the marginal log likelihood for the fixed effect parameters:
marg_nll <- function(pfixed, plist, Data, random, jnll, what = "laplace") {
  # grab the random effect vector:
  prandom <- unlist(plist[random])
  # solve the 'inner problem':
  inner <- nlminb(start = prandom, objective = jnll, pfixed = pfixed, plist = plist, Data = Data, random = random)
  # these are the empirical Bayes random effect estimates

  # calculate the Laplace approximation
  # first calculate the Hessian:
  inner$hessian <- optimHess(par = inner$par, fn = jnll, pfixed = pfixed, plist = plist, Data = Data, random = random)

  # how calculate the marginal negative log likelihood:
  inner$laplace <-
    inner$objective + # the joint log likelihood at the current set of fixed and random effects
    0.5 * log(det(inner$hessian)) - # log determinant of the Hessian
    length(prandom) / 2 * log(2 * pi) # a constant, we could ignore for optimization
  if (what == "laplace") {
    return(inner$laplace) # the marginal log likelihood
  }
  if (what == "full") {
    return(inner) # the full object so we can do things with it
  }
}

# Function to calculate finite difference gradient
# Not strictly necessary, but helpful for optimizer:
grad <- function(pfixed, ...) {
  delta <- 0.0001
  gr <- rep(0, length(pfixed))
  for (i in seq_along(gr)) {
    dvec <- rep(0, length(pfixed))
    dvec[i] <- delta
    # Calculate central finite difference
    val1 <- marg_nll(pfixed = pfixed + dvec, ...)
    val0 <- marg_nll(pfixed = pfixed - dvec, ...)
    gr[i] <- (val1 - val0) / (2 * delta)
  }
  gr
}

# Define inputs
plist <- list("logmu" = 0, "logsd" = 0, "eps" = rep(0, nrow(Data)))
random <- "eps"
fixed_pars <- unlist(plist[setdiff(names(plist), random)])

# Estimate parameters
opt2 <- nlminb(
  objective = marg_nll,
  gradient = grad,
  start = fixed_pars,
  jnll = joint_nll,
  plist = plist,
  Data = Data,
  random = random,
  control = list(trace = 1)
)
Hess <- optimHess(
  par = opt2$par,
  fn = marg_nll,
  plist = plist,
  Data = Data,
  random = random,
  jnll = joint_nll
)

final_grad <- grad(opt2$par, plist = plist, Data = Data, random = random, jnll = joint_nll)
parhat <- cbind("Estimate" = opt2$par, "SE" = sqrt(diag(solve(Hess))), "final_grad" = final_grad)
parhat

summary(fit_glmer)

# grab "empirical Bayes" random effect estimates:
out <- marg_nll(opt2$par, plist = plist, Data = Data, random = random, jnll = joint_nll, what = "full")
plot(out$par, coef(fit_glmer)$Site[,1]);abline(0, 1)
