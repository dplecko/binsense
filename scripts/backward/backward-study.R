
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

kseq <- seq.int(2L, 6L, 1L)
famseq <- c("expfam-1d")
seedseq <- seq.int(22L, 28L)
solver <- c("expfam-1d", "backward")

df <- expand.grid(kseq, famseq, seedseq, solver, stringsAsFactors = FALSE)
names(df) <- c("dimension", "family", "seed", "solver")

i_range <- seq_len(nrow(df))
evl <- list()
for (i in i_range) { 
  
  aux <- synth_data_mid(n = 10^4, k = df$dimension[i], class = df$family[i],
                        seed = df$seed[i])
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
  
  sns <- binsensate(aux$X, aux$Y, aux$R, fixy, solver = df$solver[i])
  evl[[i]] <- est_err(sns, aux)
  
  cat("\r", i)
}

vis_evl(evl, df, pz = TRUE)
