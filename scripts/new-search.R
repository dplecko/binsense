
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))

slvr <- "backward"
type <- "range"
n <- Inf
k <- 14
ngroups <- if (slvr == "em") NULL else 6
set.seed(2025)

# load data
dat <- real_data(n = n, k = k, ngroups = ngroups)

plts <- spc_lst <- list()
for (method in c("IF", "ZINF")) {
  
  for (pattern in c("agn", "x", "y")) {
    
    for (solver in slvr) {
      
      cat(method, "-", pattern, "-", solver, "\n")
      # specify the search space
      spc <- search_space(pattern = pattern, method = method, 
                          fi = seq(0, 0.2, 0.05),
                          type = type, k = k)
      
      # perform the search
      spc <- infer_search(dat, spc, solver = solver)
      spc_lst[[length(spc_lst) + 1]] <- spc
      plts[[length(plts) + 1]] <- spc_plot(spc)
    }
  }
}

# walk along a specific direction
spc_dir_if <- search_space(pattern = "x", method = "IF", 
                           fi = seq(0, 0.1, 0.025), fi2 = 0,
                           type = "range", k = k)

spc_dir_if <- infer_search(dat, spc_dir_if, solver = slvr, se = TRUE,
                           detect_viol = FALSE)

spc_dir_zinf <- search_space(pattern = "x", method = "ZINF", 
                             fi = seq(0, 0.2, 0.05), fi2 = 0,
                             type = "range", k = k)

spc_dir_zinf <- infer_search(dat, spc_dir_zinf, solver = slvr, se = TRUE,
                             detect_viol = FALSE)


# walk along the direction, with age (W) added for the EM case
if (slvr == "em") {
  
  datW <- real_data(n = n, k = k, add_W = TRUE)
  datW$W <- cbind(datW$W, datW$W^2)
  spc_dir_ifw <- search_space(pattern = "x", method = "IF", 
                              fi = seq(0, 0.1, 0.025), fi2 = 0,
                              type = "range", k = k)
  
  spc_dir_ifw <- infer_search(datW, spc_dir_ifw, solver = "em", se = TRUE,
                              detect_viol = FALSE, mc_cores = 5)
  
  spc_dir_zinfw <- search_space(pattern = "x", method = "ZINF", 
                                fi = seq(0, 0.2, 0.05), fi2 = 0,
                                type = "range", k = k)
  
  spc_dir_zinfw <- infer_search(datW, spc_dir_zinfw, solver = "em", se = TRUE,
                                detect_viol = FALSE, mc_cores = 5)
}

# check correctness with a known ground truth
spc_dir_if_gt <- infer_search(dat, spc_dir_if, solver = slvr, gt = TRUE, 
                              se = TRUE)

spc_dir_zinf_gt <- infer_search(dat, spc_dir_zinf, solver = slvr, gt = TRUE, 
                                 se = TRUE)

binsensate:::zinf_max_fi(dat$X, dat$Y, dat$R)

# save the data for plotting
# save(spc_lst, spc_dir_if, spc_dir_zinf, spc_dir_if_gt, spc_dir_zinf_gt, plts,
#      file = file.path("results", paste0(slvr, "-search.RData")))

# load(file.path("results", paste0(slvr, "-search.RData")))

# produce the plots
plts[[1]] <- spc_plot(spc_lst[[1]], spc_lst[[4]])
plts[[4]] <- NULL
plts[[6]] <- spc_plot(spc_dir_if)
cowplot::plot_grid(plotlist = plts, ncol = 3L, 
                   labels = paste0("(", letters[seq_along(plts)], ")"))

ggsave(file.path(root, "results", paste0(slvr, ".png")), width = 16, height = 7)

# ground-truth plots
ggsave(
  file.path(root, "results", paste0(slvr, "-xdir-if-gt.png")), 
  plot = spc_plot(spc_dir_if_gt),
  width = 6, height = 4
)

ggsave(
  file.path(root, "results", paste0(slvr, "-xdir-zinf-gt.png")), 
  plot = spc_plot(spc_dir_zinf_gt),
  width = 6, height = 4
)
