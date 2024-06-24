# This is an example of using the delta method to derive the standard error of
# transformation for a single parameter.

# See Jay Ver Hoef's paper:
# Who Invented the Delta Method?
# https://doi.org/10.1080/00031305.2012.687494

# Lets imagine we have some parameter named 'param'.
# We'll simulate from 'param' so we have some samples to work with for
# comparison to prove things to ourselves:

set.seed(1)
param_samples <- rnorm(1e6, 10, 0.8)
hist(param_samples, breaks = 80)

# Why a normal distribution? Because our models typically assume estimated
# parameters can be approximated well by a normal distribution. This goes
# back to the Central Limit Theorem, which essentially says the distribution
# of means from any distribution is itself normal. So, the distribution
# of a parameter estimate (which is a distribution on a mean) should
# be normally distributed given sufficient data.

# Pretend this is our parameter 'estimate' from some model:
theta <- mean(param_samples)
theta

# And pretend this is our standard error of our 'param' estimate:
theta_se <- sd(param_samples)

theta
theta_se

# Say we wanted to find the standard error of the log of our parameter.
# We will represent this transformation function as g:
# g(x) = log(x)

# We can get the mean easily:
g_theta <- log(theta)
g_theta

hist(log(param_samples), breaks = 80)
abline(v = g_theta, col = "red")

# And we can confirm this is right, because we have some samples from the parameter
# to prove it to ourselves:
mean(log(param_samples))

# But how do we find the standard error on this transformed parameter?

# We need to find the derivative of x, g'.
# This one is simple, but we could use symbolic differentiation in R:
D(expression(log(x)), "x")

# So:
# g'(x) = 1/x

# We can apply the delta method here as our derivative of the transformation
# times the original standard error:
# e.g. https://doi.org/10.1080/00031305.2012.687494 Eqn 2. (but here
# for the SD not the variance)
g_theta_se <- (1 / theta) * theta_se
g_theta_se

# Which matches our simulation-based check:
sd(log(param_samples))

# We can use that to draw our 95% CI:
.q <- qnorm(0.975)
.q
abline(v = g_theta - .q * g_theta_se, col = "blue")
abline(v = g_theta + .q * g_theta_se, col = "blue")

# The intiution (paraphrasing Wikipedia https://en.wikipedia.org/wiki/Delta_method) is that a first-order Taylor-series expansion

# Now the multivariate case (say combining multiple parameters) is a bit more
# complex because we have covariance to account for, but the principle is the
# same.
# See the file 'delta-method-multivariate.R'
