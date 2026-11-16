#set seed
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

#poission_lasso_objective
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

#test poission_lasso_objective (result aroximately 326.6079)
beta_init <- c(log(mean(y)), rep(0, n - 1))

poission_lasso_objective(beta_init, y)


#Gradient of the smoothed LASSO objective
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

#test Gradient of the smoothed LASSO objective
grad_test <- poission_lasso_grad(beta_init, y)

length(grad_test)
all(is.finite(grad_test))
grad_test[1:5]


# check the gradient using finite differences
h <- 1e-5
j <- 2

beta_plus <- beta_init
beta_minus <- beta_init

beta_plus[j] <- beta_plus[j] + h
beta_minus[j] <- beta_minus[j] - h

# Numerical gradient
grad_numeric <- (
  poission_lasso_objective(beta_plus, y) -
    poission_lasso_objective(beta_minus, y)
) / (2 * h)

# Analytical gradient
grad_analytic <- poission_lasso_grad(beta_init, y)[j]

grad_numeric
grad_analytic
abs(grad_numeric - grad_analytic)


# test the gradient at non-zero beta values
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
) / (2 * h)

grad_analytic <- poission_lasso_grad(beta_test, y)[j]

abs(grad_numeric - grad_analytic)


# Hessian of the smoothed LASSO objective
poission_lasso_hess <- function(beta, y, lambda=1, eps=1e-6) {
  
  # Calculate the linear predictor and mean
  eta <- cumsum(beta)
  mu <- exp(eta)
  n <- length(beta)
  # Initialise the Hessian matrix
  H <- matrix(0, nrow=n, ncol=n)
  # Calculate the Hessian of negative log-likelihood
  for (i in 1:n) {
    for (j in 1:n) {
      H[i,j] <- sum(mu[max(i,j):n])
    }
  }
  
  # Add the penalty second derivatives to the diagonal
  for (j in 2:n) {
    penalty_second <- lambda * eps / (beta[j]^2 + eps)^(3/2)
    H[j,j] <- H[j,j] + penalty_second
  }
  
  return(H)
}


# Check Hessian using finite differences
h <- 1e-6
j <- 2

beta_plus <- beta_test
beta_minus <- beta_test

beta_plus[j] <- beta_plus[j] + h
beta_minus[j] <- beta_minus[j] - h

# Numerical Hessian column
hess_numeric <- (
  poission_lasso_grad(beta_plus, y) -
    poission_lasso_grad(beta_minus, y)
) / (2 * h)

# Analytical Hessian column
hess_analytic <- poission_lasso_hess(beta_test, y)[,j]

# Compare the two columns
max(abs(hess_numeric - hess_analytic))

# Check all columns of the Hessian
h <- 1e-6
H_analytic <- poission_lasso_hess(beta_test, y)
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
  
  errors[j] <- max(abs(hess_numeric - H_analytic[,j]))
}

max(errors)



# fit the Poisson changepoint model using BFGS
fit <- optim(
  par = beta_init,
  fn = poission_lasso_objective,
  gr = poission_lasso_grad,
  y = y,
  lambda = 1,
  method = "BFGS",
  control = list(maxit = 1000, reltol = 1e-9)
)

# Check convergence
fit$convergence
fit$value

beta_hat <- fit$par
mu_hat <- exp(cumsum(beta_hat))

length(mu_hat)
range(mu_hat)


# plot the fitted mean function

# Plot observed data
plot(time, y, pch = 16, col = "grey",
     xlab = "Time", ylab = "Poisson count",
     main = "Poisson Changepoint Model (lambda = 1)")

# Add the true mean function
lines(time, mu_true, col = "black", lwd = 2, type = "s")
# Add the estimated mean function
lines(time, mu_hat, col = "red", lwd = 2, type = "s")
# Add a legend
legend("topright",
       legend = c("Observed data", "True mean", "Estimated mean"),
       col = c("grey", "black", "red"),
       pch = c(16, NA, NA),
       lty = c(NA, 1, 1),
       lwd = c(NA, 2, 2),
       bty = "n")
