
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

#' * phi = phi_x = phi_xy *
res <- list()
p_range <- seq(0.0, 0.3, 0.03)
df <- expand.grid(fi = p_range)
dat <- real_data()
for (j in seq_len(nrow(df))) {
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- df$fi[j]
  
  res[[j]] <- sensitivity_ATE(fixy, dat, cmb)
  df$adjusted_ATE[j] <- res[[j]][["ATE"]]
  
}

ggplot(df, aes(x = fi, y = adjusted_ATE)) +
  geom_line() + theme_bw() + geom_point() +
  xlab(latex2exp::TeX("Fidelity parameter \\phi")) +
  ylab("Adjusted ATE") +
  ggtitle(latex2exp::TeX("Basic case: $\\phi = \\phi_{x} = \\phi_{xy}$")) +
  scale_y_continuous(labels = scales::percent)

ggsave("~/Desktop/sensitivity.png", width = 6, height = 4)

#' * phi_x = phi_xy case *
df <- expand.grid(fi_x0 = p_range, fi_x1 = p_range)
df$adjusted_ATE <- 0
df$cmb_diff <- 0

for (i in seq_len(nrow(df))) {
  
  res <- sensitivity_ATE(df$fi_x1[i], df$fi_x0[i], order = 1)
  df$adjusted_ATE[i] <- res[["ATE"]]
  df$cmb_diff[i] <- res[["cmb_diff"]]
}

fidel <- matrix(df$adjusted_ATE, ncol = length(p_range))
signif <- fidel
signif[,] <- 0

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~fidel)
fig <- fig %>% add_surface(x = ~p_range, y = ~p_range, z = ~signif, opacity = 0.6)

fig <- fig %>% layout(
  scene = list(xaxis=list(title = TeX("\\phi_{x_0}")),
               yaxis=list(title = TeX("\\phi_{x_1}")),
               zaxis=list(title = TeX("ATE_{\\phi_{x_0}, \\phi_{x_1}}")))
) %>% config(mathjax = "cdn")

fig

# difference in comorbidities plot
fidel <- matrix(df$cmb_diff, ncol = length(p_range))

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(x = ~p_range, y = ~p_range, z = ~fidel)

fig <- fig %>% layout(
  scene = list(xaxis=list(title = TeX("\\phi_{x_0}")),
               yaxis=list(title = TeX("\\phi_{x_1}")),
               zaxis=list(title = "Difference in Expected # of comorbidities")
  )
)

fig # can we induce a larger difference? how? that is key for sensitivity...

# phi_xy -> need to check if it can be solved?
