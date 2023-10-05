
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

# load data
dat_lst <- real_data()
R <- dat_lst$R
X <- dat_lst$X
Y <- dat_lst$Y

# calculate correlation matrix
cor_matrix <- cor(R)

# Melt the matrix
cor_melted <- as.data.frame(as.table(cor_matrix))

# Plot
ggplot(data = cor_melted, aes(x=Var1, y=Var2)) + 
  geom_tile(aes(fill=Freq), color='white') +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
  theme_minimal() + 
  labs(fill="Correlation") + scale_fill_viridis_c() +
  geom_text(aes(label=sprintf("%.2f", Freq)), vjust=1)

des_mom0 <- r_descent_2d(var(R), R, 0, verbose = TRUE)
Sigma_hat <- tail(des_mom0, n = 1L)[[1]]

ggplot(melt(Sigma_hat), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + scale_fill_viridis_c() +
  geom_text(aes(label = round(value, 2))) +
  theme_minimal()