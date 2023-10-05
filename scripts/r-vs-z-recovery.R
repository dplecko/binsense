root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

seed_seq <- c(2022:2026)
r_err <- z_err <- rs_err <- rep(0, length(seed_seq))
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
  # des <- descent_2d(var(r), r, 0.3, n_iter = 100,
  #                   type = "Newton-Raphson", true_Sigma = attr(r, "Sigma"))
  # 
  # # vis_descent_2d(des, 10:30)
  # r_err[i] <- vec_norm(attr(r, "Sigma") - tail(des$Sigma, 1L)[[1]])
  
  # r-gradients with Method of Moments
  des_mom <- r_descent_2d(var(r), r, 0.3, n_iter = 100, multistep = TRUE)
  
  rs_err[i] <- vec_norm(attr(r, "Sigma") - tail(des_mom, 1L)[[1]])
  cat(i, "th iteration\n")
}

plot(r_err, pch = 19, ylim = c(0, max(c(r_err, rs_err), na.rm = TRUE)))
lines(r_err)
points(z_err, col="red", pch = 19)
lines(z_err, col="red")
points(rs_err, col="blue", pch = 19)
lines(rs_err, col="blue")