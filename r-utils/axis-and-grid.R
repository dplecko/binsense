
grid_to_plt <- function(lgrid, pattern, method) {
  
  is_const <- function(x) all(x == mean(x))
  
  to_model <- function(method) 
    if (method == "IM") "Independent Missingness" else "Zero-Inflated"
  
  model <- to_model(method)
  
  df <- subset(lgrid, pattern == pattern)
  txt_clr <- ifelse(df$fi_change == TRUE, "darkgray", "black")
  
  if (pattern == "agn" || 
      (pattern == "x" & (is_const(df$fi_x0y0) || is_const(df$fi_x1y0))) ||
      (pattern == "y" & (is_const(df$fi_x0y0) || is_const(df$fi_x0y1)))
     ) {
    
    if (!is.null(df$ATE_gt)) {
      
      tmp <- df$ATE
      df$ATE <- df$ATE_gt
      df$ATE_gt <- tmp
    }
    
    plt <- ggplot(df, aes(x = fi_x0y0, y = ATE, color = "Estimand")) + 
      geom_point() + geom_line() +
      theme_bw() +
      xlab(latex2exp::TeX(paste0("$\\phi$", " (", model, ")"))) + 
      ylab(latex2exp::TeX(
        "ATE Estimate"
        # "ATE$_{x_0, x_1}(y ; \\Phi = (\\phi, \\phi, \\phi, \\phi))$"
      )) +
      scale_y_continuous(labels = scales::percent)
    
    if (!is.null(df$ATE_gt)) {
      
      plt <- plt + 
        geom_line(aes(y = ATE_gt, color = "Ground Truth")) +
        geom_point(aes(y = ATE_gt, color = "Ground Truth")) +
        scale_color_manual(
          name = "Quantity",
          values = c("black", "red"), 
          labels = c("Estimated", "Ground Truth")
        )
    } else plt <- plt + scale_color_manual(values = "black") + 
                  guides(color = "none")
    
    if (!is.null(df$ATE_lwr)) {
      
      plt <- plt +
        geom_ribbon(aes(ymin = ATE_lwr, ymax = ATE_upr), alpha = 0.4)
    }
    
  } else {
    
    if (pattern == "y") df$fi_x1y0 <- df$fi_x0y1
    plt <- ggplot(df, aes(x = fi_x0y0, y = fi_x1y0)) + 
      geom_tile(aes(fill = ATE), color = 'white') +
      scale_fill_gradient2(
        low="blue", high="red", mid="white", midpoint=0,
        limits = c(-0.02, 0.02),
        labels = scales::percent
      ) +
      theme_bw() + 
      labs(fill="ATE") +
      geom_text(aes(label=sprintf("%.2f%%", 100 * ATE), color = txt_clr), vjust=1) +
      scale_color_identity() +
      xlab(latex2exp::TeX(paste0("$\\phi_{", pattern, "_0}$", " (", model, ")"))) + 
      ylab(latex2exp::TeX(paste0("$\\phi_{", pattern, "_1}$", " (", model, ")")))
  }
  
  plt +
    theme(axis.title = element_text(size = 12), 
          axis.text = element_text(size = 10))
}

grid_multi_plt <- function(grid, slvr) {

  cowplot::plot_grid(
    grid_to_plt(grid, slvr, "agnostic"), 
    grid_to_plt(grid, slvr, "x"), 
    grid_to_plt(grid, slvr, "y"),
    ncol = 3L
  )
}
