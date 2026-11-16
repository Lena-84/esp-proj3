# Practical 3: Changepoint Estimation
# Group members and contributions:
# Abdullah Alajalin (s3045135) 33%
# Xingchan Liao (s2931335) 33%
# Andreas Athinodorou (s2241333) 33%
# Github repo: https://github.com/Lena-84/esp-proj3.git

#
# Overview:
# We estimate the number and locations of changepoints in the mean of
# Poisson-distributed observations using a smoothed LASSO penalty.
# BFGS is used for model fitting, and mBIC selects the penalty parameter.
# Parametric bootstrap then assesses uncertainty, corrects estimation bias,
# and constructs 95% confidence intervals for the mean function.

############ FUNCTIONS ############

# poission_lasso_objective
# Calculate the penalised negative Poisson log-likelihood.
# beta controls the cumulative log-mean; beta[1] is not penalised.
poission_lasso_objective<- function(beta, y, lambda=1, eps=1e-6){
  #calculate the linear predictor
  eta<-cumsum(beta)
  #convert log-mean to possion mean
  mu<-exp(eta)
  #calculate nagative log-likelihood
  neg_loglik<--sum(dpois(y, lambda=mu, log=TRUE))
  #calculate the smoothed LASSO penalty
  penalty<-lambda*sum(sqrt(beta[-1]^2+eps))
  # Calculate the total penalised objective
  objective <- neg_loglik + penalty
  return(objective)
}

# Gradient of the smoothed LASSO objective
# Reverse cumulative sums give the derivative of the negative log-likelihood.
poission_lasso_grad <- function(beta, y, lambda=1, eps=1e-6) {
  # Calculate the linear predictor and mean
  eta <- cumsum(beta)
  mu <- exp(eta)
  # Gradient of negative log-likelihood
  residual <- mu - y
  grad_loglik <- rev(cumsum(rev(residual)))
  # initialize gradient of penalty
  grad_penalty <- numeric(length(beta))
  # penalty derivatives for beta[2], ..., beta[n]
  grad_penalty[-1] <- lambda * beta[-1] / sqrt(beta[-1]^2 + eps)
  # Combine both gradients
  gradient<-grad_loglik+ grad_penalty
  return(gradient)
}


# Hessian of the smoothed LASSO objective
# Each entry sums the means from max(i,j) to n.
# Reverse cumulative sums avoid repeated calculations.
poission_lasso_hess <- function(beta, y, lambda=1, eps=1e-6) {
  
  eta <- cumsum(beta)
  mu <- exp(eta)
  n <- length(beta)
  # Precompute the sum of mu from each position to n.
  remaining_sum <- rev(cumsum(rev(mu)))
  # Construct the likelihood Hessian using matrix indices.
  index <- outer(1:n, 1:n, pmax)
  H <- matrix(remaining_sum[index], nrow=n, ncol=n)
  # Add the penalty curvature to the diagonal, excluding beta[1].
  for (j in 2:n) {
    penalty_second <- lambda * eps / (beta[j]^2 + eps)^(3/2)
    H[j,j] <- H[j,j] + penalty_second
  }
  
  return(H)
}


# Use BFGS to fit the Poisson changepoint model for one penalty value.
# lambda is the penalty, y is the data and beta_init is the starting beta.
# The optim result gives the fitted coefficients in $par and the status
# in $convergence (0 means converged).
fit_lasso <- function(lambda, y, beta_init) {
  fit <- optim(
    beta_init,
    fn = poission_lasso_objective,
    gr = poission_lasso_grad,
    y = y,
    lambda = lambda,
    method = "BFGS",
    control = list(maxit = 5000, reltol = 1e-10)
  )
  return(fit)
}

############ DATA SIMULATION ############
# Simulate three Poisson-mean segments with true changes at 30 and 70.
# Set the seed so that every group member reproduces the same data.
# set seed
set.seed(12345)

n <- 100
time <- 1:n

#3 time district
tau_true <- c(1, 30, 70, n+1) 
m_true <- length(tau_true)-2
theta_true <- c(10, 2, 7)

mu_true <- numeric(n)

for(i in seq_along(theta_true))
  mu_true[tau_true[i]:(tau_true[i+1]-1)] <- theta_true[i]

y <- rpois(n, lambda = mu_true)

plot(time, y, pch = 16, col = "grey")
lines(mu_true, col = "black", lwd = 2, type = "s")
abline(v = tau_true[-c(1,length(tau_true))], lty = "dashed")

############ FUNCTION TESTS ############


# At zero jumps the fitted mean equals the sample mean.
beta_init <- c(log(mean(y)), rep(0, n-1))
cat("Initial objective:", poission_lasso_objective(beta_init, y), "\n")

# Test both zero and non-zero jump coefficients using central differences.
beta_test <- beta_init
beta_test[2] <- 0.5
beta_test[3] <- -0.3

h <- 1e-5
j <- 2
beta_plus <- beta_test
beta_minus <- beta_test
beta_plus[j] <- beta_plus[j] + h
beta_minus[j] <- beta_minus[j] - h

grad_numeric <- (
  poission_lasso_objective(beta_plus, y) -
    poission_lasso_objective(beta_minus, y)
) / (2*h)
grad_analytic <- poission_lasso_grad(beta_test, y)[j]
cat("Gradient absolute error:", abs(grad_numeric - grad_analytic), "\n")

# Verify Hessian dimensions, symmetry, finiteness, and all columns by
# central differences of the gradient.
H_test <- poission_lasso_hess(beta_test, y)
cat("Hessian dimensions:", dim(H_test), "\n")
cat("Hessian symmetric:", isSymmetric(H_test), "\n")
cat("Hessian finite:", all(is.finite(H_test)), "\n")

h <- 1e-6
errors <- numeric(n)
for (j in 1:n) {
  beta_plus <- beta_test
  beta_minus <- beta_test
  beta_plus[j] <- beta_plus[j] + h
  beta_minus[j] <- beta_minus[j] - h
  
  hess_numeric <- (
    poission_lasso_grad(beta_plus, y) -
      poission_lasso_grad(beta_minus, y)
  ) / (2*h)
  errors[j] <- max(abs(hess_numeric - H_test[,j]))
}
cat("Maximum Hessian absolute error:", max(errors), "\n")


############ INITIAL BFGS FIT ############

# Fit with lambda = 1 and check the optimiser's convergence code.
fit <- optim(
  par=beta_init,
  fn=poission_lasso_objective,
  gr=poission_lasso_grad,
  y=y,
  lambda=1,
  method="BFGS",
  control=list(maxit=1000, reltol=1e-9)
)
cat("Q3 convergence:", fit$convergence, "\n")
cat("Q3 objective:", fit$value, "\n")

beta_hat <- fit$par
mu_hat <- exp(cumsum(beta_hat))

# Compare fitted means with the true piecewise mean and observed counts.
plot(time, y, pch=16, col="grey",
     xlab="Time", ylab="Poisson count",
     main="Poisson Changepoint Model (lambda = 1)")
lines(time, mu_true, col="black", lwd=2, type="s")
lines(time, mu_hat, col="red", lwd=2, type="s")
legend("topright",
       legend=c("Observed data", "True mean", "Estimated mean"),
       col=c("grey", "black", "red"),
       pch=c(16, NA, NA),
       lty=c(NA, 1, 1),
       lwd=c(NA, 2, 2),
       bty="n")


############ PENALTY GRID SEARCH AND mBIC ############

# Start with a constant mean equal to the sample mean and zero jumps.
beta_init <- c(log(mean(y)), rep(0, n - 1))

# Grid of penalty values to try, equally spaced on the log scale.
log_lambda <- seq(-1, 5, by = 0.05)
lambda_grid <- exp(log_lambda)
n_grid <- length(lambda_grid)

# Store one result for each penalty value.
edf <- numeric(n_grid)       # effective degrees of freedom
m_tilde <- numeric(n_grid)   # effective number of changepoints
mbic <- numeric(n_grid)      # mBIC score
conv <- integer(n_grid)      # convergence code, 0 means converged

# Fit and score the model for each penalty value.
for (i in seq_len(n_grid)) {
  fit <- fit_lasso(lambda_grid[i], y, beta_init)
  conv[i] <- fit$convergence
  beta_hat <- fit$par
  
  # Evaluate both Hessians at the same fitted coefficients,
  # with the penalty (H_lambda) and without it (H_0).
  H_lambda <- poission_lasso_hess(beta_hat, y, lambda = lambda_grid[i])
  H_0 <- poission_lasso_hess(beta_hat, y, lambda = 0)
  
  # Solve H_lambda %*% A = H_0, then add the diagonal entries of A.
  # This gives the effective degrees of freedom, a measure of model
  # flexibility.
  edf[i] <- sum(diag(solve(H_lambda, H_0)))
  
  # With lambda = 0 the cost is the negative log-likelihood at beta_hat,
  # so change the sign to get the log-likelihood.
  loglik <- -poission_lasso_objective(beta_hat, y, lambda = 0)
  
  # Remove one degree of freedom for the baseline log-mean.
  m_tilde[i] <- edf[i] - 1
  
  # mBIC balances how well the model fits against its effective number
  # of changepoints.
  mbic[i] <- 3 * m_tilde[i] * log(n) - 2 * loglik
}

# optim returns convergence code 0 when a fit converges.
if (any(conv != 0)) {
  warning(sum(conv != 0), " grid fits did not converge")
}

# The best lambda is the one with the smallest mBIC.
i_best <- which.min(mbic)
lambda_hat <- lambda_grid[i_best]


############ SELECTED MODEL AND CHANGEPOINTS ############

# Refit the model using the selected penalty.
fit_best <- fit_lasso(lambda_hat, y, beta_init)
if (fit_best$convergence != 0) {
  warning("best fit did not converge")
}
beta_best <- fit_best$par

# The running sum gives the log-mean at each time point.
mu_hat <- exp(cumsum(beta_best))

# Compare the fitted mean with the true mean and the observed counts.
plot(
  time, y, pch = 16, col = "grey", xlab = "Time", ylab = "y",
  main = "True and estimated mean function"
)
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_hat, col = "red", lwd = 2, type = "s")
legend(
  "topright",
  legend = c("Data", "True mean", "Estimated mean"),
  col = c("grey", "black", "red"),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 1),
  lwd = c(NA, 2, 2)
)

# Round down the effective number of changepoints.
m_hat <- floor(m_tilde[i_best])

# Skip the baseline and rank the jumps by absolute size. Add 1 to recover
# the original time indices, then sort the selected locations.
jump_order <- order(abs(beta_best[-1]), decreasing = TRUE) + 1
tau_hat <- sort(jump_order[seq_len(m_hat)])

# Mark the estimated changepoints on the plot.
abline(v = tau_hat, col = "red", lty = 2)

cat("Best lambda =", round(lambda_hat, 3),
    "(log lambda =", log_lambda[i_best], ")\n")
cat("Estimated number of changepoints =", m_hat, "\n")
cat("Estimated changepoint locations =", tau_hat, "\n")



############ PARAMETRIC BOOTSTRAP AND CONFIDENCE INTERVALS ############