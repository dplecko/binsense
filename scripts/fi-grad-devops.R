
k <- 5
pz <- runif(2^k)
pz <- pz / sum(pz)
pos <- sample(2^k, 1)

infer_pz(pos, pz, 0.3, k) - infer_pz(pos, pz, rep(0.3, k), k)

# get an example that goes out of the box
dat_lst <- real_data()
k <- 5
R <- dat_lst$R[, seq_len(k)]
X <- dat_lst$X
Y <- dat_lst$Y

backward_solver(X, Y, R, list(list(0.2, 0.2), list(0.2, 0.2)))
