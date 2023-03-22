
seed_seq <- c(2022:2031)
r_err <- z_err <- rep(0, length(seed_seq))
k <- 6
for (i in seq_along(seed_seq)) {
  
  # z-gradients
  z <- z_2d(k, 10^4, seed = seed_seq[i])
  des <- descent_2d(var(z), z, 0, n_iter = 100,
                    type = "Newton-Raphson", true_Sigma = attr(z, "Sigma"))
  
  # vis_descent_2d(des, 1:100)
  z_err[i] <- vec_norm(attr(z, "Sigma") - tail(des$Sigma, 1L)[[1]])
  
  # r-gradients
  r <- r_2d(k, 10^4, 0.3, seed = seed_seq[i])
  des <- descent_2d(var(r), r, 0.3, n_iter = 100,
                    type = "Newton-Raphson", true_Sigma = attr(r, "Sigma"))
  
  # vis_descent_2d(des, 10:30)
  r_err[i] <- vec_norm(attr(r, "Sigma") - tail(des$Sigma, 1L)[[1]])
  cat(i, "th iteration\n")
}

plot(r_err, pch = 19, ylim = c(0, max(r_err, na.rm = TRUE)))
lines(r_err)
points(z_err, col="red", pch = 19)
lines(z_err, col="red")

setwd("~/ICU/op")
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))


### Hessian computation
k <- 3
Sigma <- gen_Sigma(k) # array(0, dim = c(k, k)) #
pz <- cpp_scaling_2d(Sigma)

z <- z_2d(3, 10^4, Sigma = Sigma)
des <- descent_2d(var(z), z, 0, n_iter = 10,
                  type = "Newton-Raphson", true_Sigma = attr(z, "Sigma"))

vis_descent_2d(des, 1:10)
hess <- cpp_hessian_2d(pz / sum(pz), which(!upper.tri(Sigma))-1, k)

solve(hess, rep(1, 6))

### check if gradients match! ()
set.seed(2023)
k <- 3
z <- z_2d(k, 10^4, seed = 2023)
Sigma <- attr(z, "Sigma")

gradient_2d_fi0(Sigma, z)
gradient_2d(Sigma, z, fi = 0) #' * tick *

### check if likelihoods match
likelihood_2d(Sigma, z, fi = 0)
likelihood_2d_fi0(Sigma, z)

# test gradient descent with and w/o Hessian

t(z) %*% Sigma
dim(t(z))

### perform gradient descent
k <- 3
fi <- 0.1
r <- r_2d(k, 10^4, fi, seed = 2028)
des <- descent_2d(var(r), r, fi, n_iter = 10,
                  type = "Newton-Raphson", 
                  true_Sigma = attr(r, "Sigma"))

vis_descent_2d(des, 1:100)
z_err[i] <- vec_norm(attr(z, "Sigma") - tail(des$Sigma, 1L)[[1]])


#' * Order of Terms: * -> gradient exploration
seed_seq <- c(2022:2031)
kseq <- c(2, 3, 4, 5, 6)
fiseq <- c(0, 0.01, 0.025, 0.05, 0.1, 0.2)
res <- expand.grid(seed = seed_seq, k = kseq, fi = fiseq)
res$gam_t <- res$gam_m <- 0
res_lst <- list()

for (i in seq_len(nrow(res))) {
  
  fi <- res$fi[i]
  seed <- res$seed[i]
  k <- res$k[i]
  
  # r-gradients
  r <- r_2d(k, 10^4, fi, seed = seed)
  des <- descent_2d(var(r), r, fi, n_iter = 20,
                    type = "Newton-Raphson", true_Sigma = attr(r, "Sigma"))
  
  res_lst[[i]] <- des
  # vis_descent_2d(des, 10:30)
  res$gam_t[i] <- mean(tail(des$gsc))
  res$gam_m[i] <- max(des$gsc)
  cat(i, "th iteration\n")
}

ggplot(reshape2::melt(res, id.vars = names(res)[1:3]), 
       aes(x = fi, y = value, color = variable)) +
  geom_point() + geom_smooth() +
  facet_grid(cols = vars(k)) + theme_bw()
