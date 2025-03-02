
test_that("ill posedness helpers", {
  
  set.seed(2024)
  
  dat <- synth_data_mid(n = 10^2, k = 5,
                        fi = list(list(0.01, 0.01), list(0.01, 0.01)),
                        Sigma = matrix(runif(25, min = -1, max = 0.5), ncol = 5),
                        class = "expfam-2d")
  
  add_rows <- 20
  dat$R <- rbind(dat$R, matrix(1, nrow = add_rows, ncol = 5))
  dat$X <- c(dat$X, rep(1, add_rows))
  dat$Y <- c(dat$Y, rep(1, add_rows))
  
  suppressMessages(
    expect_message(
      binsensate(dat$X, dat$Y, dat$R, fi = 0.3, method = "IF", 
                 A = NULL, solver = "backward"),
      regexp = "fi parameter was updated"
    )
  )
})
