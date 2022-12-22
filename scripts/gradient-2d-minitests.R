
# gen_Sigma works
Sigma <- gen_Sigma(2)

# scaling_2d works
scaling_2d(Sigma)

# z_2d works
table(z_2d(2, 10^4) %*% c(1, 2)) / sum(table(z_2d(2, 10^4) %*% c(1, 2)))

attr(z_2d(2, 10^4), "pz") / sum(attr(z_2d(2, 10^4), "pz"))


# r_2d works
table(r_2d(2, 10^4, 0.3) %*% c(1, 2)) / sum(table(r_2d(2, 10^4, 0.3) %*% c(1, 2)))

attr(r_2d(2, 10^4, 0.3), "pz") / sum(attr(r_2d(2, 10^4, 0.3), "pz")) * (1 - 0.3)^2


# get_idx works (tick)
get_idx(c(0, 1), 1, 1)

# get_fi_hash works
get_fi_hash(c(0, 0), 0.3)

# compute_grad_ij works
fi_wgh_hash <- list()

compute_grad_ij(1, 2, c(4), c(0.25, 0.25, 0.25, 0.25), 0.3)
