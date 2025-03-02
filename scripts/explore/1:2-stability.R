

A <- matrix(
  c(
    1/2+eps, 0, 1/2-eps,
    0, 1/2+eps, 1/2-eps,
    1/2-eps, 0, 1/2+eps
  ), ncol = 3, byrow = TRUE
)

eigen(A)

v <- c(2, 1, -5/2) * 2 / 3 / sqrt(5)
as.numeric(t(v) %*% A %*% v)

library(progress)


v <- c(1, 2, -2)
v <- v / vec_norm(v)
t(v) %*% A %*% v

nrep <- 100000
n <- 2
eps <- 0.01
pb <- progress_bar$new(
  format = "(:spin) [:bar] :percent",
  total = nrep, clear = FALSE, width = 60
)

for (i in 1:nrep) {
  
  A <- matrix(runif(n * n), nrow = n)
  A <- A - diag(diag(A))
  
  A <- A / (2 * rowSums(A)) * (1-eps)
  A <- A + diag(rep((1+eps) / 2, n))
  
  if (any(Re(eigen(A)$values) <= 0)) browser()
  pb$tick()
}

# sup vs. eigenvalue
stoch_mat <- function(n) {
  mat <- matrix(runif(n * n), nrow = n)
  mat <- mat / rowSums(mat)
  return(mat)
}


nrep <- 100
d <- 8
gamma <- c()
for (i in seq_len(nrep)) {
  
  A <- stoch_mat(d)
  Ainv <- solve(A)
  gamma <- c(gamma, max(abs(Ainv)) / max(abs(eigen(Ainv)$value)))
}

max(gamma)
