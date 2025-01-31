
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))

library(parallel)
library(progress)

n_seq <- c(1000, 2000) # c(1000, 5000, 10000)
k_seq <- c(3, 5) # c(3, 5, 10)
nseed <- 50
ncores <- parallel::detectCores() / 2

method <- "ZINF"
ret <- c()

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) elapsed: :elapsed",
  total = length(n_seq) * length(k_seq) * nseed, clear = FALSE, width = 60
)

for (n in n_seq) {
  
  for (k in k_seq) {
    
    for (seed in seq_len(nseed)) {
      
      set.seed(seed)
      nparams <- if (method == "IM") k else 2^k
      fi <- list(
        list(runif(nparams, min = 0, max = 1/3), runif(nparams, min = 0, max = 1/3)), 
        list(runif(nparams, min = 0, max = 1/3), runif(nparams, min = 0, max = 1/3))
      )
      
      # generate the data
      dat <- synth_data_mid(n = n, k = k, seed = seed, fi = fi, 
                            lam = runif(k, min = -2/3, max = 2/3),
                            mu = runif(k, min = -2/3, max = 2/3),
                            beta = runif(1, min = -0.4, max = 0),
                            icept_x = runif(1, min = -0.5, max = 0.5),
                            icept_y = runif(1, min = -0.5, max = 0.5),
                            class = "expfam-2d")
      
      # get the Louis standard error
      beta_new <- binsensate(X = dat$X, Y = dat$Y, R = dat$R, fi = fi, 
                             solver = "two-stage-em", method = method,
                             new = TRUE)$beta
      beta_old <- binsensate(X = dat$X, Y = dat$Y, R = dat$R, fi = fi, 
                             solver = "two-stage-em", method = method, 
                             new = FALSE)$beta
      
      ret <- rbind(ret, data.frame(n = n, k = k, seed = seed,
                                   sign = beta_new > beta_old, 
                                   diff = beta_new - beta_old,
                                   beta_new = beta_new, beta_old = beta_old))
      
      pb$tick()
      
      # cat("Finished iteration n =", n, "; k =", k, "; seed =", seed, "\n")
    }
  }
}

ret <- as.data.table(ret)

ggplot(
  data = ret
) +
  geom_density(aes(x = beta_new - beta_old, fill = "black"), alpha = 0.5) +
  # geom_density(aes(x = beta_new, fill = "black"), alpha = 0.5) +
  # geom_density(aes(x = beta_old, fill = "red"), alpha = 0.5) +
  facet_grid(rows = vars(n), cols = vars(k))

ret[, median(beta_new - beta_old), by = c("n", "k")]
