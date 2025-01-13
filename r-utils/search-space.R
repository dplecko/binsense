
rng_to_df <- function(pattern, fi_seq) {
  
  pattern <- match.arg(pattern, c("agn", "x", "y"))
  if (pattern == "agn") {
    
    df <- data.frame(pattern = rep(pattern, length(fi_seq)))
    df$fi_x0y0 <- df$fi_x1y0 <- df$fi_x0y1 <- df$fi_x1y1 <- fi_seq
  } else if (pattern == "x") {
    
    fi_grid <- expand.grid(a = fi_seq, b = fi_seq)
    
    df <- data.frame(pattern = rep(pattern, nrow(fi_grid)))
    df$fi_x0y0 <- df$fi_x0y1 <- fi_grid$a
    df$fi_x1y0 <- df$fi_x1y1 <- fi_grid$b
  } else if (pattern == "y") {
    
    fi_grid <- expand.grid(a = fi_seq, b = fi_seq)
    
    df <- data.frame(pattern = rep(pattern, nrow(fi_grid)))
    df$fi_x0y0 <- df$fi_x1y0 <- fi_grid$a
    df$fi_x0y1 <- df$fi_x1y1 <- fi_grid$b
  }
  
  df
}

search_space <- function(pattern = c("agn", "x", "y"),
                         type = c("range", "bound"),
                         method = c("IM", "ZINF"),
                         fi, k) {
  
  pattern <- match.arg(pattern, c("agn", "x", "y"))
  type <- match.arg(type, c("range", "bound"))
  method <- match.arg(method, c("IM", "ZINF"))
  
  params <- rng_to_df(pattern, fi_seq = fi)
  params$ATE <- NA
  
  if (method == "ZINF") k <- 2^k
  
  if (type == "bound") {
    
    opt_params <- lapply(
      1:nrow(params),
      function(i) {
        
        fi <- list(
          list(rep(params[i,]$fi_x0y0, k), rep(params[i,]$fi_x0y1, k)),
          list(rep(params[i,]$fi_x1y0, k), rep(params[i,]$fi_x1y1, k))
        )
      }
    )
  } else opt_params <- NULL
  
  list(
    pattern = pattern, fi = fi, type = type, method = method,
    params = params, opt_params = opt_params
  )
}

infer_search <- function(dat, spc, solver, nbreaks = 5) {
  
  pattern <- spc$pattern
  if (spc$type == "bound") {
    
    for (i in seq_len(nrow(spc$params))) {
      
      flag <- TRUE
      while (flag) {
        
        flag <- FALSE
        for (x in c(0, 1)) {
          
          if ((pattern == "agn" || pattern == "y") & x == 1) next
          for (y in c(0, 1)) {
            
            if ((pattern == "agn" || pattern == "x") & y == 1) next
            for (j in seq_along(spc$opt_params[[i]][[1+x]][[1+y]])) {
              
              # get the current pattern
              fi_curr <- spc$opt_params[[i]]
              
              # get the (x, y) subpattern and its j-th value
              fij_curr <- fi_curr[[1+x]][[1+y]][j]
              
              # define the search range and number of breakpoints
              fi_range <- seq(0, spc$params[i,][[paste0("fi_x", x, "y", y)]], 
                              length.out = nbreaks)
              
              if (length(fi_range) == 0) next
              
              # prepare for search
              fij_opt <- NA
              f_opt <- -Inf
              for (l in seq_len(nbreaks)) {
                
                cand_fij <- fi_range[l]
                
                # evaluate the candidate
                fi_test <- fi_curr
                if (pattern == "agn") {
                  
                  for (xp in c(0, 1)) for (yp in c(0, 1))
                    fi_test[[1+xp]][[1+yp]][j] <- cand_fij
                } else if (pattern == "x") {
                  
                  for (yp in c(0, 1))
                    fi_test[[1+x]][[1+yp]][j] <- cand_fij
                } else if (pattern == "y") {
                  
                  for (xp in c(0, 1))
                    fi_test[[1+xp]][[1+y]][j] <- cand_fij
                }
                
                f_tst <- binsensate(dat$X, dat$Y, dat$R, fi = fi_test, 
                                    method = spc$method, solver = solver)$ATE
                
                if (f_tst > f_opt) {
                  
                  
                  f_opt <- f_tst
                  fij_opt <- cand_fij
                  spc$opt_params[[i]] <- fi_test
                }
              }
              
              # if an update is found, renew the flag
              if (fij_opt != fij_curr) {
                
                flag <- TRUE
                cat("Found a pattern improvement.",
                    "fi_old =", fij_curr, "vs.", "fi_new =", fij_opt, "\n")
              }
            }
          }
        }
      }
    }
  }
  
  for (i in seq_len(nrow(spc$params))) {
    
    if (spc$type == "bound") {
      
      fi <- spc$opt_params[[i]]
    } else {
      
      fi <- list(
        list(spc$params[i,]$fi_x0y0, spc$params[i,]$fi_x0y1),
        list(spc$params[i,]$fi_x1y0, spc$params[i,]$fi_x1y1)
      )
    }
    
    res <- binsensate(dat$X, dat$Y, dat$R, fi = fi, method = spc$method,
                      solver = solver)
    spc$params$ATE[i] <- res$ATE
  }
  
  spc
}

spc_plot <- function(spc) {
  
  pattern <- spc$pattern
  grid_to_plt(spc$params, spc$pattern, spc$method)
}
