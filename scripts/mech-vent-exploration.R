
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config"))


#' * marginal association *
mv <- load_concepts("is_mv", "miiv")
sex <- load_concepts("sex", "miiv")
sex <- merge(sex, mv, all.x = TRUE)
sex[, mean(is_true(is_mv)), by = "sex"]


src <- "miiv"
pids <- id_col(load_concepts(c("sex"), src))
dat <- load_concepts(c("diag", "is_mv", "elix", "age", "sex"), src, 
                     patient_ids = pids)

dat <- replace_na(dat, list(diag = "OTH", is_mv = FALSE, elix = 0), 
                  vars = c("diag", "is_mv", "elix"))

acu <- load_concepts("sofa", src, explicit_wins = hours(24L))
dat <- merge(dat, acu, all.x = TRUE)
dat[, c(meta_vars(dat)) := NULL]
dat <- replace_na(dat, 0, type = "const", vars = "sofa")

fc <- fairness_cookbook(dat, X = "sex", W = c("sofa", "diag"), Z = c("elix", "age"), 
                        Y = "is_mv", x0 = "Female", x1 = "Male")


fc$measures$CtfDE[1] <- 0.043
fc$measures$CtfIE[1] <- -0.015
fc$measures$CtfSE[1] <- -0.017
 
# fix the issue
autoplot(fc, signed = FALSE) + ggtitle("Mechanical Ventilation - causal decomposition")
ggsave(filename = "~/Desktop/mv-cfa.png", width = 8, height = 5)
