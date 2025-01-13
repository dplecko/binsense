
test_that("binsensate handles different input types", {
  
  set.seed(2024)
  
  # create the A matrix
  A <- matrix(0, nrow = 2^3, ncol = 2^3)
  A[upper.tri(A)] <- runif(sum(upper.tri(A)), max = 0.04)
  diag(A) <- 1 - rowSums(A)
  
  tst_set <- list(
    A1 = list(k = 3, fi = 0.25, method = "IM", A = NULL),
    A2 = list(k = 3, fi = c(0.2, 0.1, 0.05), method = "IM", A = NULL),
    A3 = list(
      k = 3,
      fi = list(
        list(c(0.1, 0.1, 0.05), c(0.2, 0.1, 0.05)),
        list(c(0.03, 0.1, 0.05), c(0.2, 0.0, 0.1))
      ),
      method = "IM", A = NULL
    ),
    B1 = list(k = 3, fi = 0.1, method = "ZINF", A = NULL),
    B2 = list(k = 3, fi = runif(2^3, 0.1), method = "ZINF", A = NULL),
    B3 = list(
      k = 3,
      fi = list(
        list(runif(2^3, 0.05), runif(2^3, 0.05)),
        list(runif(2^3, 0.05), runif(2^3, 0.05))
      ),
      method = "ZINF", A = NULL
    ),
    C = list(k = 3, fi = NULL, method = NULL, A = A)
  )
  
  for (i in seq_along(tst_set)) {
    
    st <- tst_set[[i]]
    dat <- synth_data_mid(n = 10^3, k = 3,
                          fi = st$fi, method = st$method, A = st$A,
                          class = "expfam-2d")
    
    for (solver in c("backward", "two-stage-em")) {
      result <- binsensate(dat$X, dat$Y, dat$R, fi = st$fi, method = st$method, 
                           A = st$A, solver = solver)
      expect_true(is.numeric(result$ATE))
    }
  }
})
