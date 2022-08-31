
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

sens_manual <- function(dat, fi_x1, fi_x0, cmb) {
  
  add_cmb <- function(col, bmi, fi_x1, fi_x0) {

    add <- rbinom(n = length(col), size = 1, 
                  prob = ifelse(bmi < 25, fi_x0, fi_x1))
    col | add
  }
  
  for (col in cmb) {
    
    dat[, c(col) := add_cmb(get(col), get("bmi"), fi_x1, fi_x0)]
  }
  
  logreg <- glm(
    as.formula(paste("factor(death) ~", paste(c("(bmi >25)", cmb), collapse = "+"))),
    data = dat, family = "binomial"
  )
  
  cmb_diff <- dat[, list(n_cmb = sum(.SD)), by = c("stay_id", "bmi"), .SDcols = cmb]
  diff <- tapply(cmb_diff$n_cmb, cmb_diff$bmi < 25, mean)

  list(
    ATE = mean(predict(logreg, dat[, bmi := 30], type = "response")) - 
          mean(predict(logreg, dat[, bmi := 20], type = "response")),
    cmb_diff = diff[2] - diff[1]
  )
} 


p_range <- seq(0, 0.3, 0.02)
df <- expand.grid(fi_x0 = p_range, fi_x1 = p_range)
df$adjusted_ATE <- 0
df$cmb_diff <- 0

for (i in seq_len(nrow(df))) {
  
  cat("\r", round(i / nrow(df) * 100), "%")
  res <- sens_manual(copy(dat), df$fi_x1[i], df$fi_x0[i], cmb)
  df$adjusted_ATE[i] <- res[["ATE"]]
  df$cmb_diff[i] <- res[["cmb_diff"]]
}

fidel <- matrix(df$adjusted_ATE, ncol = length(p_range))
signif <- fidel
signif[,] <- 0

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(x = ~p_range, y =~p_range, z = ~fidel)
fig <- fig %>% add_surface(x = ~p_range, y =~p_range, z = ~signif, opacity = 0.6)

fig <- fig %>% layout(
  scene = list(xaxis=list(title = TeX("AX")),
               yaxis=list(title = TeX("\\phi_{x_1}")),
               zaxis=list(title = TeX("ATE_{\\phi_{x_0}, \\phi_{x_1}}")))
) %>% config(mathjax = "cdn")

fig

# difference in comorbidities plot
fidel <- matrix(df$cmb_diff, ncol = length(p_range))

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(x = ~p_range, y =~p_range, z = ~fidel)

fig <- fig %>% layout(
  scene = list(xaxis=list(title = "Fidelity x0"),
  yaxis=list(title = "Fidelity x1"),
  zaxis=list(title = "Difference in Expected # of comorbidities")
  )
)

fig # can we induce a larger difference? how? that is key for sensitivity...




