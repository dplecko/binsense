
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
    l1_rel = abs(sns$ATE - aux$ATE) / abs(aux$ATE),
    pz_MSE = sum( (sns$pz - aux$pz)^2 ),
    runtime = sns$runtime / 10^9
  )
}

vis_evl <- function(evl, df, type = c("ate", "pz", "comp", "l1_rel", "all"), 
                    split = TRUE, out_df = FALSE) {
  
  type <- match.arg(type, c("ate", "pz", "comp", "l1_rel", "all"))
  
  df <- cbind(df, ATE_MSE = do.call(c, lapply(evl, `[[`, "ATE_MSE")))
  df <- cbind(df, pz_MSE = do.call(c, lapply(evl, `[[`, "pz_MSE")))
  df <- cbind(df, runtime = do.call(c, lapply(evl, `[[`, "runtime")))
  df <- cbind(df, l1_rel = do.call(c, lapply(evl, `[[`, "l1_rel")))
  
  if (out_df) return(df)
  
  p1 <- ggplot(df, aes(x = dimension, y = ATE_MSE, color = solver, fill = solver)) +
    geom_point() + geom_smooth(method = "loess", formula = "y ~ x") + theme_bw() +
    xlab("Dimensionality") + ylab("MSE(ATE)") +
    theme(legend.position = "bottom")
  
  p2 <- ggplot(df, aes(x = dimension, y = pz_MSE, color = solver, fill = solver)) +
    geom_point() + geom_smooth(method = "loess", formula = "y ~ x") + theme_bw() +
    xlab("Dimensionality") + ylab("MSE(P(z))") +
    theme(legend.position = "bottom")
  
  p3 <- ggplot(df, aes(x = dimension, y = runtime, color = solver, fill = solver)) +
    geom_point() + geom_smooth(method = "loess", formula = "y ~ x") + theme_bw() +
    xlab("Dimensionality") + ylab("binsensate() runtime") +
    theme(legend.position = "bottom")
  
  p4 <- ggplot(df, aes(x = dimension, y = l1_rel, color = solver, fill = solver)) +
    geom_point() + geom_smooth(method = "loess", formula = "y ~ x") + theme_bw() +
    xlab("Dimensionality") + ylab("Relative L1 error") +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "bottom")
  
  if (split) {
    
    p1 <- p1 + facet_grid(cols = vars(family))
  }
  
  switch(type,
         ate = p1,
         pz = p2,
         comp = p3,
         l1_rel = p4,
         all = cowplot::plot_grid(p1, p2, p3, p4, ncol = 2L))

}

