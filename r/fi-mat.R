
fi_to_A_xy <- function(fi, k, method) {
  
  if (k < 10) {
    
    switch (method,
            IF = cpp_A_xy_if(fi, k),
            ZINF = cpp_A_xy_zinf(fi, k)
    )
  } else { # use sparse matrices for large k
    
    switch (method,
            IF = cpp_A_xy_if_sparse(fi, k),
            ZINF = cpp_A_xy_zinf_sparse(fi, k)
    )
  }
}

fi_to_Ainv_xy <- function(fi, k, method) {
  
  switch (method,
          IF = cpp_Ainv_xy_if(fi, k),
          ZINF = cpp_Ainv_xy_zinf(fi, k)
  )
}

fi_to_A <- function(fi, k, method) {
  
  A <- list(list(NULL, NULL), list(NULL, NULL))
  for (x in c(0, 1)) for (y in c(0, 1)) {
    
    A[[1+x]][[1+y]] <- fi_to_A_xy(fi[[1+x]][[1+y]], k, method)
  }
  
  A
}

expand_fi <- function(fi, k, method) {
  
  if (method == "ZINF") k <- 2^k
  
  # expand fi if needed
  for (x in c(T, F))
    for (y in c(T, F))
      if (length(fi[[x + 1]][[y + 1]]) == 1 & k > 1) 
        fi[[x + 1]][[y + 1]] <- rep(fi[[x + 1]][[y + 1]], k)
  
  fi
}

expit <- function(x) exp(x) / (1 + exp(x))

vec_norm <- function(x) sqrt(sum(x^2))