
# taskset -c 0-63 Rscript scripts/explore/boot-vs-louis.R
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))

library(parallel)
library(progress)

n_seq <- c(1000, 2000, 5000) # c(1000, 5000, 10000)
k_seq <- c(3, 5, 7) # c(3, 5, 10)
nboot <- 64
nseed <- 50

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) elapsed: :elapsed",
  total = length(n_seq) * length(k_seq) * nseed, clear = FALSE, width = 60
)

bvl <- c()
for (n in n_seq) {
  
  for (k in k_seq) {
    
    for (seed in seq_len(nseed)) {
      
      set.seed(seed)
      fi <- list(
        list(runif(k, min = 0, max = 1/3), runif(k, min = 0, max = 1/3)), 
        list(runif(k, min = 0, max = 1/3), runif(k, min = 0, max = 1/3))
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
      frun <- binsensate(X = dat$X, Y = dat$Y, R = dat$R, fi = fi, 
                         solver = "em", se = TRUE)
      
      se_louis <- frun$beta_se
      
      # get the boostrap estimate of the error
      beta_bvals <- mclapply(
        seq_len(nboot), 
        function(boot) {
          b_idx <- sample.int(n = n, replace = TRUE)
          brun <- binsensate(X = dat$X[b_idx], Y = dat$Y[b_idx], R = dat$R[b_idx, ], 
                             fi = fi, solver = "em")
          brun$beta
        }, mc.cores = n_cores()
      )
      
      se_boot <- sd(do.call(c, beta_bvals))
      bvl <- rbind(bvl, data.frame(n = n, k = k, seed = seed,
                                   beta = frun$beta, se_louis = se_louis,
                                   se_boot = se_boot))
      
      pb$tick()
    }
  }
}

# visualization
bvl <- as.data.table(bvl)
bvl[, LBR := se_louis / se_boot] # LBR < 1 indicates undercoverage of Louis (too narrow)

ggplot(bvl, aes(x = LBR)) +
  geom_density() +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", linewidth = 1) +
  facet_grid(rows = vars(n), cols = vars(k)) +
  theme_bw()
