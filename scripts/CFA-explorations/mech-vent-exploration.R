
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), source))
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
dat <- replace_na(dat, 0, type = "const", vars = "sofa")


# Y (outcome) disparity analysis
dat1 <- data.table::copy(dat)
dat1 <- dat1[, c(meta_vars(dat1)) := NULL]
fc <- fairness_cookbook(dat1, X = "sex", W = c("sofa", "diag"), Z = c("elix", "age"), 
                        Y = "is_mv", x0 = "Female", x1 = "Male")


# fc$measures$CtfDE[1] <- 0.043
# fc$measures$CtfIE[1] <- -0.015
# fc$measures$CtfSE[1] <- -0.017
 
# fix the issue
autoplot(fc, signed = FALSE) + ggtitle("Mechanical Ventilation - causal decomposition")
ggsave(filename = "~/Desktop/mv-cfa.png", width = 8, height = 5)

#' * compute the \Delta on this example *

# add death
dod <- load_concepts("death", src)
dat <- merge(dat, dod, all.x = TRUE)
dat[is.na(death), death := FALSE]

cdat <- data.table::copy(dat)

cdat[, sex := as.integer(sex == "Female")]
cdat <- cbind(cdat, model.matrix(~., data = cdat[, "diag"])[, -1])
cdat[, diag := NULL]

diag_cols <- grep("diag", names(cdat), value = TRUE)
crf <- causal_forest(
  W = cdat$is_mv, # treatment
  Y = cdat$death, # outcome
  X = cdat[, c(diag_cols, "elix", "age", "sex", "sofa"), with = FALSE] # covariates
)

# very small values; check marginal differences
dat$Delta <- crf$predictions # need to just check the sign

fc_Delta <- fairness_cookbook(dat, X = "sex", W = c("sofa", "diag"), 
                              Z = c("elix", "age"), Y = "Delta", x0 = "Female", 
                              x1 = "Male")

autoplot(fc_Delta) + ggtitle(TeX("$\\Delta$-Mech. Vent. decomposition on MIMIC-IV"))

