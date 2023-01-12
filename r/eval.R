
# evl_frm <- function(fwd, aux) {
#   
#   data.frame(adist = unlist(fwd$adverse), 
#              bdist = abs(unlist(fwd$beta) - aux$beta),
#              zdist = unlist(
#                lapply(fwd$pz, function(z) 100 * sum(abs(z - aux$pz)))
#              ),
#              iter = seq_along(fwd$adverse))
# }

est_err <- function(sns, aux) {
  
  list(
    ATE_MSE = min(abs(sns$ATE - aux$ATE)^2, 1),
    pz_MSE = sum( (sns$pz - aux$pz)^2 )
  )
}

vis_evl <- function(evl, df, pz = FALSE) {
  
  df <- cbind(df, ATE_MSE = do.call(c, lapply(evl, `[[`, "ATE_MSE")))
  df <- cbind(df, pz_MSE = do.call(c, lapply(evl, `[[`, "pz_MSE")))
  
  p1 <- ggplot(df, aes(x = dimension, y = ATE_MSE, color = solver, fill = solver)) +
    geom_point() + geom_smooth(method = "loess", formula = "y ~ x") + theme_bw() +
    xlab("Dimensionality") + ylab("MSE(ATE)")
  
  p2 <- ggplot(df, aes(x = dimension, y = pz_MSE, color = solver, fill = solver)) +
    geom_point() + geom_smooth(method = "loess", formula = "y ~ x") + theme_bw() +
    xlab("Dimensionality") + ylab("MSE(P(z))")
  
  if (pz) {
    
    cowplot::plot_grid(p1, p2, ncol = 1L)
  } else {
    
    p1
  }

}

