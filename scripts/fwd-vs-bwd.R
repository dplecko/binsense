

k <- 2
fi <- 0.5
type <- "rand"

if (type == "equal") {
  pz <- rep(2^(-k), 2^k)
} else {
  breaks <- runif(2^k - 1)
  breaks <- c(0, sort(breaks), 1)
  pz <- diff(breaks)
}

A <- array(0, dim = c(2^k, 2^k))

for (i in seq_len(2^k)) {
  
  for (j in seq_len(2^k)) {
    
    i_bin <- as.integer(intToBits(i-1))
    j_bin <- as.integer(intToBits(j-1))
    
    if (any(i_bin > j_bin)) next
    if (all(i_bin == j_bin)) next
    # if (sum(j_bin - i_bin) > order) next
    A[i, j] <- fi^sum(j_bin - i_bin)
  }
}
A <- A + diag(1 - colSums(A))

pr <- A %*% pz

# perform marginal inference
pri <- rep(0, k)
for (i in seq_len(k)) {
  
  pri[i] <- 0
  
  for (j in seq_len(2^k)) {
    
    if (as.integer(intToBits(j-1))[i] == 1) pri[i] <- pri[i] + pr[j]
    
  }
  
}

pzi <- pri / (1 - fi)

trzi <- (1 - pzi) / (1-pzi + fi * pzi)

B <- array(0, dim = c(2^k, 2^k))
for (i in seq_len(2^k)) {
  
  for (j in seq_len(2^k)) {
    
    i_bin <- as.integer(intToBits(i - 1))
    j_bin <- as.integer(intToBits(j - 1))
    
    if (any(i_bin < j_bin)) next
    if (all(i_bin == j_bin)) next
    
    idx <- ((i_bin - j_bin) > 0)[seq_len(k)]
    
    B[i, j] <- prod(1 - trzi[idx]) * prod(trzi[!idx])
  }
}

B <- B + diag(1 - colSums(B))

cat("L2 norm between A^(-1) and B", sum((solve(A) - B)^2), "\n")

cat("L2 norm between Bp(r) and p(z)", sum((B %*% pr - pz)^2), "\n")

cat("L2 norm between A^(-1)p(r) and p(z)", sum((solve(A) %*% pr - pz)^2), "\n")

cat("L2 norm between A^(-1)p(r) and Bp(r)", sum((solve(A) %*% pr - B %*% pr)^2), 
    "\n")

