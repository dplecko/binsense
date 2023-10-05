
library(reticulate)

# Specify the path to the Python script
python_script <- file.path(root, "basicpy.py")

# Source the Python script
source_python(python_script)

# generate some R
k <- 9
fi <- 0.1 
r <- r_2d(k, 10^4, fi, seed = 2028)

fi_mat <- array((1-fi)^2, dim = c(k, k))
diag(fi_mat) <- 1-fi
mu_hat <- t(r) %*% r / nrow(r) / fi_mat

Sigma_zero <- array(0, dim = c(k, k))

# perform descent with R
microbenchmark::microbenchmark(
  Sigma_hat <- mom_gradient(Sigma_zero, mu_hat, -1, nsamp = 500), 
  times = 1L
)

# which is better?
microbenchmark::microbenchmark(
  Sigma_hat <- mom_gradient(Sigma_zero, mu_hat, -1, nsamp = 500, T, F), 
  times = 1L
)

Sigma_hat <- mom_gradient(Sigma_zero, mu_hat, -1, nsamp = 500)

cat("L2-norm of R implementation", sum((attr(r, "Sigma") - Sigma_hat)^2))

Sigma_pyhat <- mom_grad_py(mu_hat, num_iterations = 500L)

cat("L2-norm of R implementation", sum((attr(r, "Sigma") - Sigma_pyhat)^2))
