
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))

### ground truth testing
n <- Inf
k <- 14
datW <- real_data(n = n, k = k, add_W = TRUE)

fi <- list(list(0.05, 0.05), list(0, 0))
bns <- binsensate(datW$X, datW$Y, datW$R, fi, W = datW$W, method = "IF",
                  se = TRUE, mc_cores = 5, tau_xw = function(w) w^2)

# spc_dir_ifw <- search_space(pattern = "x", method = "IF", 
#                             fi = seq(0, 0.1, 0.025), fi2 = 0,
#                             type = "range", k = k)
# 
# spc_dir_ifw <- infer_search(datW, spc_dir_ifw, solver = "em", se = TRUE,
#                             detect_viol = FALSE, mc_cores = 5,
#                             tau_xw = function(w) w^2)
# save(spc_dir_ifw, file = file.path(root, "age_real_only.RData"))

load(file.path(root, "age_real_only.RData"))

p_w <- spc_plot(spc_dir_ifw)

spc_dir_ifw_gt <- infer_search(datW, spc_dir_ifw, solver = "em", gt = TRUE, 
                               se = TRUE, mc_cores = 5, tau_xw = function(w) w^2)

p_wgt <- spc_plot(spc_dir_ifw_gt)

### age plot
data <- real_data(add_W = TRUE, ret_dt = TRUE)
data[, bmi_bin := (bmi_all >= 25)]

library(mgcv)
m <- gam(bmi_bin ~ s(age), family = binomial, 
         data = data[, c("bmi_bin", "age"), with=FALSE])
pse <- predict(m, newdata = data.frame(age = 18:91), se.fit = TRUE)
gam_vals <- data.frame(
  px = expit(pse$fit), px_lwr = expit(pse$fit - 1.96 * pse$se.fit), 
  px_upr = expit(pse$fit + 1.96 * pse$se.fit), age = 18:91
)

raw_vals <- data[, .(px = mean(bmi_bin), 
                     se = sqrt(mean(bmi_bin) * (1 - mean(bmi_bin)) / .N)), 
                 by = age]

p_age <- ggplot(gam_vals, aes(x = age, y = px)) +
  geom_line(color = "orange", linewidth = 1.2) +
  geom_ribbon(aes(ymin = px_lwr, ymax = px_upr), alpha = 0.1) +
  theme_bw() + xlab("W (age in years)") + ylab("P(X = 1 | W = w)") +
  geom_point(data=raw_vals) +
  scale_y_continuous(labels = scales::percent)

save(p_age, p_w, p_wgt, spc_dir_ifw, spc_dir_ifw_gt,
     file = file.path(root, "age_final.RData"))

# load("~/ICU/op/age_final.RData")

cowplot::plot_grid(p_age, p_w, labels = c("(a)", "(b)"))

ggsave("results/age_analysis.png", width = 10, height = 4)
