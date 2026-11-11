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