
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))

devtools::load_all()

dat_lst <- real_data(add_W = TRUE)
k <- 5
nboot <- 50
samp <- TRUE
R <- dat_lst$R[samp, seq_len(k)]
X <- dat_lst$X[samp]
Y <- dat_lst$Y[samp]
W <- dat_lst$W[samp]

fi <- 0.1
if (!is.list(fi)) fi <- list(list(fi, fi), list(fi, fi))
fi <- expand_fi(fi, ncol(R), "IF")
A <- fi_to_A(fi, ncol(R), "IF")

solved <- em_solver_cts(X, Y, R, as.matrix(W), fi, "IF", A, se = TRUE)
solved_d <- em_solver(X, Y, R, fi, "IF", A, se = TRUE)

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

load("~/ICU/op/age_fisearch.RData")
p_incw <- spc_plot(spc_dir_ifw)

cowplot::plot_grid(p_age, p_incw, labels = c("(a)", "(b)"))

ggsave("results/age_analysis.png", width = 10, height = 4)

### ground truth testing
n <- 10^5
k <- 3
datW <- real_data(n = n, k = k, add_W = TRUE)
datW$W <- cbind(datW$W, datW$W^2)

spc_dir_ifw <- search_space(pattern = "x", method = "IF", 
                            fi = seq(0, 0.1, 0.025), fi2 = 0,
                            type = "range", k = k)

spc_dir_ifw <- infer_search(datW, spc_dir_ifw, solver = "em", se = TRUE,
                            detect_viol = FALSE, mc_cores = 5)

spc_dir_ifw_gt <- infer_search(datW, spc_dir_ifw, solver = "em", gt = TRUE, 
                               se = TRUE)

spc_plot(spc_dir_ifw_gt)

# extract a case where there is a clear mismatch
fi = list(list(0.2, 0.2), list(0, 0))
n = 10^5
k = 3
d = 2
method <- "IF"
lambda <- c(-0.156770960233408, -0.398341578744453, 0.134068566391964, 
            -0.0462034360129478, -0.26682004624256)
icept_x <- 0.8380211

mu <- c(1.62570429560644, 0.79891976300749, 0.993580973877462, 0.525550538227234, 
        -0.0506855754074046)
icept_y <- -4.282786 

Sigma <- structure(c(-1.41566621505028, 0.586336536513105, 0.403890292748151, 
                     0.586336536513105, -3.07200575746431, -0.0496813156739413, 0.403890292748151, 
                     -0.0496813156739413, -1.75957570947528), dim = c(3L, 3L))
SLO <- structure(c(-1.41566621505028, 0.586336536513105, 0.403890292748151, 
                   0.260796362451022, 0.366268025912868, 0.586336536513105, -3.07200575746431, 
                   -0.0496813156739413, 0.192504785514286, 0.0647035367408105, 0.403890292748151, 
                   -0.0496813156739413, -1.75957570947528, 0.860289045427417, 0.630535963143475, 
                   0.260796362451022, 0.192504785514286, 0.860289045427417, 1.25387745028644, 
                   0.310919340180732, 0.366268025912868, 0.0647035367408105, 0.630535963143475, 
                   0.310919340180732, 0.651545090323566), dim = c(5L, 5L))

beta <- 0.2775879

# intervene in parameters to make Z _||_ W
# SLO[4:5, 1:3] <- SLO[1:3, 4:5] <- 0

pz <- binsensate:::theta_Ez(SLO, 4:5) / sum(binsensate:::theta_Ez(SLO, 4:5))
colSums(pz * (binsensate:::cpp_idx_to_bit(1:2^k - 1, dim = k) %*% 
                t(solve(SLO[4:5, 4:5]) %*% SLO[4:5, 1:3])))


seed <- 2022
dati <- synth_data_mid(
  n = n, k = k,
  seed = seed, fi = fi, method = method,
  class = "expfam-2d-cts",
  Sigma = Sigma, 
  SLO = SLO,
  lam = lambda, 
  mu = mu,
  icept_x = icept_x, icept_y = icept_y,
  beta = beta
)

bns <- binsensate(dati$X, dati$Y, dati$R, fi, W = dati$W, method = method)
SLO - tail(bns$SLO, n = 1)[[1]]
tail(bns$LMB, n = 1)[[1]]

# is the data correct w/o measurement error?
SLO - tail(binsensate:::cormat_to_Sigma_cts(t(cbind(dati$Z, dati$W)) %*% cbind(dati$Z, dati$W) / nrow(dati$Z), k, d), n = 1L)[[1]]


mod <- glm(
  Y ~ ., data = data.frame(Y = dati$Y, X = dati$X, Z = dati$Z, W = dati$W), family = "binomial"
)

bns_now <- binsensate(dati$X, dati$Y, dati$R, fi, method = method, solver = "em")


# is W drifting?
dati_2 <- synth_data_mid(
  n = n, k = k,
  seed = seed, fi = fi, method = method,
  class = "expfam-2d-cts",
  Sigma = NULL, 
  SLO = tail(bns$SLO, n = 1)[[1]],
  lam = c(tail(bns$LMB, n = 1)[[1]]$lambda_z, tail(bns$LMB, n = 1)[[1]]$lambda_w), 
  mu = c(tail(bns$LMB, n = 1)[[1]]$mu_z, tail(bns$LMB, n = 1)[[1]]$mu_w),
  icept_x = tail(bns$LMB, n = 1)[[1]]$lambda_icept, 
  icept_y = tail(bns$LMB, n = 1)[[1]]$mu_icept,
  beta = beta
)
