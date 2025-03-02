
library(umap)
library(ggplot2)
library(ggtext)
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))
options(scipen = 999)

# load the real data
dat_lst <- real_data()
k <- 11
samp <- TRUE
R <- dat_lst$R[samp, seq_len(k)]
X <- dat_lst$X[samp]
Y <- dat_lst$Y[samp]
fixy <- list(list(0.1, 0.1), list(0.1, 0.1))
init <- list()

for (i in 1:50) {
  
  res <- binsensate(X, Y, R, fi = fixy, solver = "em", n_epoch = 10,
                    rand_init = TRUE)
  
  init[[i]] <- list(
    res$theta[[1]], res$theta[[length(res$theta)]]
  )
}

data <- do.call(rbind, lapply(seq_along(init), function(i) {
  t0 <- t(do.call(c, init[[i]][[1]]))  # Convert first list to a single row
  t10 <- t(do.call(c, init[[i]][[2]])) # Convert second list to a single row
  colnames(t0) <- colnames(t10) <- NULL
  rbind(data.frame(id=paste0("$T_", i, "^(0)$"), t0),
        data.frame(id=paste0("$T_", i, "^(10)$"), t10))
}))

data$group <- as.integer(grepl("10", data$id))

res <- NULL
rm_cols <- which(colnames(data) %in% c("id", "group"))
for (i in seq_len(nrow(data))) {
  
  for (j in seq(i+1, nrow(data))) {
    
    dist <- vec_norm(
      unlist(data[i, -rm_cols]) - unlist(data[j, -rm_cols])
    )
    
    res <- rbind(
      res,
      c(dist, data$group[i], data$group[i] == data$group[j])
    )
    
  }
  
}

res <- as.data.frame(res)
names(res) <- c("dist", "group", "within")

ggplot(res, aes(x = dist, fill = interaction(group, within))) +
  geom_density() + theme_bw()

# Preparing data
mds_cols <- paste0("X", seq_along(unlist(init[[1]][[1]])))
data_mds <- data[, mds_cols]

# MDS
dist_matrix <- dist(data_mds)
mds_result <- cmdscale(dist_matrix, k = 2)

# Convert to data frame for plotting
mds_df <- as.data.frame(mds_result)
names(mds_df) <- c("Dim1", "Dim2")
mds_df$group <- data$group  # Add group for coloring

ggplot(mds_df, aes(x = Dim1, y = Dim2, fill = as.factor(group))) +
  geom_point(shape = 21, color = "black", size = 2.5) +
  theme_bw() +
  labs(x = "Dimension 1", y = "Dimension 2") +
  scale_fill_discrete(
    name = "Parameter Group",
    labels = c(
      latex2exp::TeX("$\\theta_i^{(0)}$ (Initial)"),
      latex2exp::TeX("$\\theta_i^{(10)}$ (Converged)")
    )
  ) +
  theme(legend.position = "bottom", legend.box.background = element_rect(),
        legend.text = element_text(size = 12))

ggsave(file.path(root, "paper", "figures", "init-sense.png"),
       width = 6, height = 4)
