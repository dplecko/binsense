
# test A_xy

fi <- 0.3

for (k in 3:10) {
  
  if (max(abs(A_xy(fi, k) - cpp_A_xy(fi, k))) > 10^(-10)) {
    cat("Found a mistake\n")
  }
}


A_xy(fi, k)
cpp_A_xy(fi, k)

# test invert_A_xy

for (k in 3:10) {
  
  A <- A_xy(fi, k)
  pz <- runif(2^k)
  pr <- runif(2^k)
  
  if (max(abs(invert_A_xy(A, pz, pr) - cpp_invert_A_xy(A, pz, pr))) > 10^(-10)) {
    cat("Found a mistake\n")
  }
  
  
}

invert_A_xy(A, pz, pr)
cpp_invert_A_xy(A, pz, pr)

k <- 4

for (i in seq_len(2^k)) {
  
  
  if (any(get_idx(to_bits(i-1, k), 1, 1) != cpp_get_idx(i-1, 1, 1, k))) {
    cat("Mistake for i =", i, "\n")
    break
  }
  
}


for (i in seq_len(2^k)) {
  
  print(
    get_fi_hash(to_bits(i-1, k), 0.3) - cpp_fi_hash(i-1, 0.3, k)
  )
  
}

k <- 7

microbenchmark::microbenchmark(
  lapply(seq_len(2^k),
         function(x) get_fi_hash(to_bits(i-1, k), 0.3)
  ),
  times = 5L
)

microbenchmark::microbenchmark(
  lapply(seq_len(2^k),
         function(x) cpp_fi_hash(i-1, 0.3, k)
  ),
  times = 5L
)


setwd("~/ICU/op")
source("scripts/forward/second-order-gradient.R")
k <- 7
fi <- 0.3
n_samp <- 10000
r <- r_2d(k, n_samp, fi, seed = 23)
rhash <- 1 + rowSums(r * t(replicate(nrow(r), 2^(seq_len(k) - 1))))

# cpp_grad_ij(rhash, 1, 2, k, 0.3, rep(1, 2^k))

microbenchmark::microbenchmark(
  lapply(seq_len(20),
         function(x) compute_grad_ij(1, 2, r, rep(1, 2^k), 0.3, 1)
  ),
  times = 5L
)

microbenchmark::microbenchmark(
  lapply(seq_len(20),
         function(x) compute_grad_ij(1, 2, r, rep(1, 2^k), 0.3, 6)
  ),
  times = 5L
)

max(abs(
  compute_grad_ij(1, 2, r, rep(1, 2^k), 0.3) - 
    compute_grad_ij_old(1, 2, r, rep(1, 2^k), 0.3)
))
