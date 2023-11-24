
binsensate_boot <- function(
    X, Y, R, fi, 
    solver = c("backward", "expfam-1d", "expfam-2d",
               "expfam-1d-direct", "backward-direct",
               "expfam-2d-mom", "expfam-2d-mom-grad",
               "two-stage-em"),
    mc_samples = 10^4, mc_per_samp = 10, 
    ret_em = FALSE, nboot = 100, ...
  ) {
  
  solver <- match.arg(solver, c("backward", "expfam-1d", "expfam-2d", 
                                "expfam-1d-direct", "backward-direct",
                                "expfam-2d-mom", "expfam-2d-mom-grad",
                                "two-stage-em"))
  
  bt_obj <- list()
  for (i in seq_len(nboot)) {
   
    if (i == 1) idx <- seq_along(X) else {
      
      idx <- sample(length(X), replace = TRUE)
    }
    
    Xb <- X[idx]
    Yb <- Y[idx]
    Rb <- R[idx, ]
    
    bt_obj[[i]] <- switch(
      solver,
      backward = backward_solver(Xb, Yb, Rb, fi, ...),
      `expfam-1d` = expfam_1d_solver(Xb, Yb, Rb, fi, mc_samples),
      `expfam-2d` = expfam_2d_solver(Xb, Yb, Rb, fi, mc_samples, ret_em = ret_em),
      `expfam-1d-direct` = expfam_1d_direct(Xb, Yb, Rb, fi, mc_samples),
      `expfam-2d-mom` = expfam_2d_mom(Xb, Yb, Rb, fi, mc_samples),
      `expfam-2d-mom-grad` = expfam_2d_mom_grad(Xb, Yb, Rb, fi, mc_samples),
      `two-stage-em` = two_stage_em(Xb, Yb, Rb, fi, mc_per_samp)
    ) 
    
  }
  
  bt_obj
}