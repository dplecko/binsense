
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))

n <- 10000 # 121119
k <- 5
type <- "range"

# load data
dat <- real_data(n = n, k = k)
plts <- list()
for (method in c("IM", "ZINF")) {
  
  for (pattern in c("agn", "x", "y")) {
    
    for (solver in c("two-stage-em")) {
      
      cat(method, "-", pattern, "-", solver, "\n")
      # specify the search space
      spc <- search_space(pattern = pattern, method = method, 
                          fi = seq(0, 0.2, 0.05),
                          type = type, k = k)
      
      # perform the search
      spc <- infer_search(dat, spc, solver = solver)
      plts[[length(plts) + 1]] <- spc_plot(spc)
    }
  }
}

cowplot::plot_grid(plotlist = plts, ncol = 3L)

