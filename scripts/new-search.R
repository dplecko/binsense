
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))

slvr <- "backward"
type <- "range"
n <- 121119
k <- if (slvr == "two-stage-em") 5 else 5

# load data
dat <- real_data(n = n, k = k)

binsensate(dat$X, dat$Y, dat$R, fi = 0.25, method = "ZINF",
           solver = "two-stage-em", detect_viol = FALSE)

plts <- list()
for (method in c("IM", "ZINF")) {
  
  for (pattern in c("agn", "x", "y")) {
    
    for (solver in slvr) {
      
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

# walk along a specific direction, and check the ground truth
spc_dir <- search_space(pattern = "x", method = "IM", 
                    fi = seq(0, 0.2, 0.05), fi2 = 0,
                    type = "range", k = 5)

spc_dir <- infer_search(dat, spc_dir, solver = slvr)
spc_plot(spc_dir)

spc_dir <- infer_search(dat, spc_dir, solver = slvr, gt = TRUE)
spc_plot(spc_dir)
