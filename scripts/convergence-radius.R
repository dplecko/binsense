
n <- 5
A <- matrix(rnorm(n^2), ncol = n)
x <- rnorm(n)

solve(A, x)

eigen(A, only.values = TRUE)$values

des$lam
des$mu

k <- 4
fi <- 0.05
r <- r_2d(k, 10^4, fi, seed = seed_seq[i])
des <- descent_2d(var(r), r, fi, n_iter = 20,
                  type = "Newton-Raphson", true_Sigma = attr(r, "Sigma"))
vis_conv_radius(des)
