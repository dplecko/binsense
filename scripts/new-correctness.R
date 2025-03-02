
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))

set.seed(2024)

# create the A matrix
A <- matrix(0, nrow = 2^3, ncol = 2^3)
A[upper.tri(A)] <- runif(sum(upper.tri(A)), max = 0.04)
diag(A) <- 1 - rowSums(A)

tst_set <- list(
  A1 = list(k = 3, fi = 0.25, method = "IF", A = NULL),
  A2 = list(k = 3, fi = c(0.2, 0.1, 0.05), method = "IF", A = NULL),
  A3 = list(
    k = 3,
    fi = list(
      list(c(0.1, 0.1, 0.05), c(0.2, 0.1, 0.05)),
      list(c(0.03, 0.1, 0.05), c(0.2, 0.0, 0.1))
    ),
    method = "IF", A = NULL
  ),
  B1 = list(k = 3, fi = 0.1, method = "ZINF", A = NULL),
  B2 = list(k = 3, fi = runif(2^3, 0.1), method = "ZINF", A = NULL),
  B3 = list(
    k = 3,
    fi = list(
      list(runif(2^3, max = 0.1), runif(2^3, max = 0.1)),
      list(runif(2^3, max = 0.1), runif(2^3, max = 0.1))
    ),
    method = "ZINF", A = NULL
  ),
  C = list(k = 3, fi = NULL, method = NULL, A = A)
)

nrep <- 10
res <- c()
for (i in seq_along(tst_set)) {
  
  st <- tst_set[[i]]
  
  for (j in seq_len(nrep)) {
    
    dat <- synth_data_mid(n = 10^3, k = 3, seed = j,
                          fi = st$fi, method = st$method, A = st$A,
                          class = "expfam-2d")
    
    for (solver in c("backward", "em")) {
      result <- binsensate(dat$X, dat$Y, dat$R, fi = st$fi, method = st$method, 
                           A = st$A, solver = solver)
      
      if (result$ATE > 1 || result$ATE < -1) browser()
      res <- rbind(res, data.frame(ate = result$ATE, ate_true = dat$ATE,
                                   scm = i, rep = j, solver = solver))
    }
  }
}

res <- as.data.table(res)
res[, ate_true := median(ate_true), by = "scm"]
ggplot(res, aes(x = solver, y = ate, fill = solver)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = ate_true), color = "red", linetype = "dashed") +
  facet_wrap(~ scm, scales = "free") +
  labs(
    x = "Solver",
    y = "ATE Estimate",
    fill = "Solver",
    title = "Boxplot of ATE Estimates by Solver and SCM"
  ) +
  theme_minimal()
