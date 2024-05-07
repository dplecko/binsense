
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

kseq <- seq.int(5L, 8L, 1L)
famseq <- c("expfam-2d")
seedseq <- seq.int(35L, 39L)
solver <- c("backward", "expfam-2d-mom-grad", "two-stage-em")

df <- expand.grid(kseq, famseq, seedseq, solver, stringsAsFactors = FALSE)
names(df) <- c("dimension", "family", "seed", "solver")

i_range <- seq_len(nrow(df))
evl <- list()
for (i in i_range) { 
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.1
  
  aux <- synth_data_mid(n = 10000, k = df$dimension[i], class = df$family[i],
                        fi = fixy, seed = df$seed[i])
  
  if (df$solver[i] == "two-stage-em") mc_samp <- 10 else mc_samp <- 10^4
  mcb <- microbenchmark::microbenchmark(
    sns <- binsensate(aux$X, aux$Y, aux$R, fixy, solver = df$solver[i],
                      mc_samples = mc_samp), 
    times = 1L
  )
  sns$runtime <- mcb$time
  
  evl[[i]] <- est_err(sns, aux)
  
  cat("\r", i)
}

vis_evl(evl, df, type = "all")
