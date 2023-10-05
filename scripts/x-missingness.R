
#' * Complete cases vs. Consider BMI missingness *
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

# what kind of insight can be obtained?

# use the typical case

# compare new and old Elixhauser scores
data <- load_concepts(c("death", "bmi_all", "elix_wide", "age"), "miiv", 
                      id_type = "patient")
data <- replace_na(data, val = FALSE, vars = "death")
data[, c(index_var(data)) := NULL]
data[, Obesity := NULL]
cmb_sel <- c("CHF", "Arrhythmia", "PVD", "Paralysis", "NeuroOther", "Pulmonary", 
             "Liver", "Lymphoma", "Mets", "Coagulopathy", "FluidsLytes")

data <- replace_na(data, val = 0, vars = cmb_sel)
data <- data[, c("death", "bmi_all", "age", cmb_sel), with=FALSE]


# MCAR - complete case analysis
cdat <- data[complete.cases(data)]

Zc <- as.matrix(cdat[, cmb_sel, with=FALSE])
Xc <- as.integer(cdat$bmi_all > 25)
Yc <- as.integer(cdat$death)

cmod <- glm(Y ~ ., data = data.frame(Z = Zc, X = Xc, Y = Yc), 
            family = "binomial")
summary(cmod)

# what is the ATE when you consider complete cases?
mean(
  predict(cmod, newdata = data.frame(Z = Zc, X = 1, Y = Yc), type = "response") -
    predict(cmod, newdata = data.frame(Z = Zc, X = 0, Y = Yc), type = "response")
)

# MAR - full data analysis (but without imputation)
Z <- as.matrix(data[, cmb_sel, with=FALSE])
X <- as.integer(data$bmi_all > 25)
Y <- as.integer(data$death)

# get P(y | z)
pyz <- glm(Y ~ Z, family = "binomial")

# get P(x* | y, z, R_x = 0) on complete data
p_xsyz <- glm(X ~ ., data = data.frame(Z = Zc, X = Xc, Y = Yc), 
              family = "binomial")

# predict on full data
pxs_y1 <- predict(p_xsyz, newdata = data.frame(Z = Z, Y = 1), type = "response")
pxs_y0 <- predict(p_xsyz, newdata = data.frame(Z = Z, Y = 0), type = "response")

# compute the P(y | x, z)
py_x1z <- pxs_y1 * pyz$fitted.values / (
  pxs_y0 * (1 - pyz$fitted.values) + pxs_y1 * pyz$fitted.values
)

py_x0z <- (1 - pxs_y1) * pyz$fitted.values / (
  (1 - pxs_y0) * (1 - pyz$fitted.values) + (1 - pxs_y1) * pyz$fitted.values
)

mean((py_x1z / (1 - py_x1z)) / py_x0z / (1 - py_x0z))
