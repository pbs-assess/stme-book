# Here we will use the delta method to calculate a standard error
# on a derived quantity involving 2 parameters:

# To make our code more readable:
invlogit <- plogis

# Imagine we have a bunch of fish lengths:
x <- seq(1, 20, length.out = 400)

# And a true intercept and slope parameter from which we will simulate:
b0 <- -3
b1 <- 0.4

# Form the linear predictor:
mu <- invlogit(b0 + b1 * x)
plot(x, mu, ylim = c(0, 1))

# Apply observation error:
set.seed(1)
y <- rbinom(length(x), size = 1, prob = mu)

# Plot it:
plot(x, y)

# Fit a GLM to it:
fit <- glm(y ~ x, family = binomial())
coef(fit)
est <- predict(fit) # linear predictor

# What's the length at 50% maturity?
# We can calculate that as:
calc_length_p50 <- function(x1, x2, prob_threshold = 0.5) {
  -(log((1/prob_threshold) - 1) + x1) / x2
}
(b0_hat <- coef(fit)[[1]])
(b1_hat <- coef(fit)[[2]])
p50 <- calc_length_p50(b0_hat, b1_hat)
p50

# Plot it:
plot(x, invlogit(est), ylim = c(0, 1), type = "l")
abline(v = p50, col = "red")
abline(h = 0.5, lty = 2)

# OK, but what about the standard error on p50?
# We can use the multivariate delta method.

# Grab our covariance matrix:
cov <- vcov(fit)

# Define our transformation of interest:
g <- ~ -(log((1/0.5) - 1) + b0_hat) / b1_hat

# Create a function that calculates the partial derivatives:
deriv_function <- deriv(g, c("b0_hat", "b1_hat"))
deriv_function

# Evaluate that function:
e <- eval(deriv_function)

# Note that this gives us the same thing as we did above to
# calculate p50, but it also gives us the partial derivatives:
e

# We can grab those partial derivatives:
J <- attr(e, "gradient")

# These are also known as the Jacobian: a matrix of partial derivatives of the
# entries in g(theta) with respect to the entries in theta itself (theta is our
# parameter vector).

# And matrix multiply the Jacobian with the original covariance matrix by the
# transposed Jacobian to get the variance of the derived quantity:
variance <- J  %*% cov %*% t(J )
variance

# The standard error is the square root of that:
se <- sqrt(diag(new_cov))
se

# Plot it:
abline(v = p50 + 1.96 * se, col = "blue")
abline(v = p50 - 1.96 * se, col = "blue")

