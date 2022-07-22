
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

# Question 1 - how many unique R's even appear?
ind <- load_concepts(c("bmi", "death"), "miiv")
ind[, c(index_var(ind)) := NULL]
ind <- ind[!is.na(bmi)]

dat <- merge(ind, get_icd(index = "charlson"), all.x = TRUE)
dat[is.na(death), death := FALSE]

keep_vars <- c("stay_id", "bmi", "death")
for (col in setdiff(names(dat), keep_vars)) {
  dat[is.na(get(col)), c(col) := FALSE]
}
# cmb <- setdiff(names(dat), keep_vars)
cmb <- c("DM", "Cancer", "LiverSevere", "Renal", "Stroke", "MI")
enc_mat <- t(replicate(nrow(dat), 2^(seq_along(cmb) - 1L)))

# Question 2 - can we get the A matrix?
k <- length(cmb)

# Question: are we walking out of the simplex?
# Question 3 - can we fit a model for mortality?
logreg <- glm(
  as.formula(paste("factor(death) ~", paste(c("(bmi >25)", cmb), collapse = "+"))),
  data = dat, family = "binomial"
)

# summary(logreg)
pyrx0 <- pyrx1 <- rep(0, 2^k)
logreg$coefficients[seq_len(k)+2] <- 3
for (i in seq_len(2^k)) {
  
  i_bin <- as.integer(intToBits(i-1))
  pyrx0[i] <- expit(sum(logreg$coefficients * c(1, 0, i_bin)))
  pyrx1[i] <- expit(sum(logreg$coefficients * c(1, 1, i_bin)))
  
}

prx0 <- tabulate(rowSums(dat[bmi<= 25, c(cmb), with=F] * enc_mat) + 1)
prx0 <- c(prx0, rep(0L, 2^k - length(prx0))) # zero-padding
prx0 <- prx0 / sum(prx0)

prx1 <- tabulate(rowSums(dat[bmi > 25, c(cmb), with=F] * enc_mat) + 1)
prx1 <- c(prx1, rep(0L, 2^k - length(prx1))) # zero-padding
prx1 <- prx1 / sum(prx1)

sensitivity_ATE <- function(fi_x1, fi_x0, order = Inf) {
  
  Ax0 <- Ax1 <- array(0, dim = c(2^k, 2^k))
  
  for (i in seq_len(2^k)) {
    
    for (j in seq_len(2^k)) {
      
      i_bin <- as.integer(intToBits(i-1))
      j_bin <- as.integer(intToBits(j-1))
      
      if (any(i_bin > j_bin)) next
      if (all(i_bin == j_bin)) next
      if (sum(j_bin - i_bin) > order) next
      Ax0[i, j] <- fi_x0^sum(j_bin - i_bin)
      Ax1[i, j] <- fi_x1^sum(j_bin - i_bin) 
    }
  }
  
  Ax0 <- Ax0 + diag(1 - colSums(Ax0))
  Ax1 <- Ax1 + diag(1 - colSums(Ax1))
  
  if (any(Ax0 < 0) | any(Ax1 < 0)) return(NA)
  
  if (is.element(0, eigen(Ax0)$values)) browser()
  pzx0 <- solve(Ax0) %*% prx0
  if (check_simplex(pzx0)) {

    pzx0 <- locate_simplex(prx0, solve(Ax0))
  }
  pzx1 <- solve(Ax1) %*% prx1
  if (check_simplex(pzx1)) {

    pzx1 <- locate_simplex(prx1, solve(Ax1))
    # pzx1 <- search_lambda(prx1, solve(Ax1))
  }
  if (any(c(range(pzx0), range(pzx1)) > 1) | 
      any(c(range(pzx0), range(pzx1)) < 0)) {
    cat("Walking out of simplex", c(range(pzx0), range(pzx1)), "\n")
  }
  pz <- solve(Ax0) %*% prx0 * mean(dat$bmi <= 25) + 
    solve(Ax1) %*% prx1 * mean(dat$bmi > 25) 
  
  dlt_x0x1 <- rep(0, 2^k)
  for (i in seq_len(2^k)) {
    
    dlt_x0x1[i] <- sum(pyrx1 * Ax1[, i] - pyrx0 * Ax0[, i])
  }
  
  sum(dlt_x0x1 * pz)
}

p_range <- seq(0, 0.2, 0.01)
df <- expand.grid(fi_x0 = p_range, fi_x1 = p_range)
df$adjusted_ATE <- 0

for (i in seq_len(nrow(df))) {
  
  df$adjusted_ATE[i] <- sensitivity_ATE(df$fi_x1[i], df$fi_x0[i], order = 1)
}

fidel <- matrix(df$adjusted_ATE, ncol = length(p_range))
signif <- fidel
signif[,] <- 0

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~fidel)
fig <- fig %>% add_surface(z = ~signif, opacity = 0.6)

fig <- fig %>% layout(
  scene = list(xaxis=list(title = "Fidelity x0"),
               yaxis=list(title = "Fidelity x1"),
               zaxis=list(title = "Treatment Effect Estimate",
                          range = c(-0.061, 0.01)))
)

fig
