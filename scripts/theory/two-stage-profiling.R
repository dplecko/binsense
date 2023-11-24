
library(microbenchmark)
library(profvis)

root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

# load data
dat_lst <- real_data()
R <- dat_lst$R
X <- dat_lst$X
Y <- dat_lst$Y

fixy <- list(list(0.1, 0.1), list(0.1, 0.1))

microbenchmark(
  binsensate(X, Y, R, fixy, solver = "two-stage-em", nepoch = 3),
  times = 1L
)

# the cpp version is to be updated
microbenchmark(
  binsensate(X, Y, R, fixy, solver = "two-stage-em", mc_twister = "cpp"),
  times = 1L
)

profvis({
  # Your R code here
  two_stage_em(X, Y, R, fixy, 10, n_epoch = 1)
})

Rprof("myprofile.out")

two_stage_em(X, Y, R, fixy, 10, n_epoch = 1)

Rprof(NULL)

summaryRprof("myprofile.out")

# conclusion: 1 iteration of EM takes about 1 minute
# speeding up likely not possible: profiling suggests most of the time is
# already spent in the C++ code