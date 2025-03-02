

test_that("binsensate handles chunking correctly", {
  
  set.seed(2024)
  
  dat <- synth_data_mid(n = 10^2, k = 6, fi = 0.1, method = "IF",
                        class = "expfam-2d")
  
  add_rows <- 500
  dat$R <- rbind(dat$R, matrix(0, nrow = add_rows, ncol = 6))
  dat$X <- c(dat$X, rep(1, add_rows))
  dat$Y <- c(dat$Y, rep(1, add_rows))
  
  expect_message(
    binsensate(dat$X, dat$Y, dat$R, fi = 0.1, method = "IF", 
               solver = "em"),
    regexp = "Hessian Singular"
  )
})
