
test_that("fi_mat sparse works", {
  library(Matrix)
  expect_true(
    inherits(binsensate:::fi_to_A_xy(rep(0.1, 11), 11, "IF"), "dgCMatrix")
  )
  expect_true(
    inherits(binsensate:::fi_to_A_xy(rep(0.01, 2^11), 11, "ZINF"), "dgCMatrix")
  )
})
