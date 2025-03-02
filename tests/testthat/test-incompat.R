
test_that("incompatibility of ", {
  
  set.seed(2024)
  
  for (k in c(3, 6)) {
    
    dat <- synth_data_mid(n = 10^3, k = k,
                          fi = list(list(0.01, 0.01), list(0.01, 0.01)),
                          Sigma = matrix(runif(k^2, min = -1, max = 0.5), ncol = k),
                          class = "expfam-2d")
    
    add_rows <- 20
    dat$R <- rbind(dat$R, matrix(1, nrow = add_rows, ncol = k))
    dat$X <- c(dat$X, rep(1, add_rows))
    dat$Y <- c(dat$Y, rep(1, add_rows))
    
    suppressMessages(
      expect_message(
        binsensate(dat$X, dat$Y, dat$R, fi = 0.3, method = "ZINF", 
                   A = NULL, solver = "em", detect_viol = TRUE),
        regexp = "problem may be ill-posed"
      )
    )
    
    suppressMessages(
      expect_message(
        binsensate(dat$X, dat$Y, dat$R, fi = 0.3, method = "ZINF", 
                   A = NULL, solver = "backward"),
        regexp = "fi parameter was updated"
      )
    )
  }
})
