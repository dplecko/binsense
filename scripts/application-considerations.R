
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

# for MIMIC-IV, OMR data is loaded on id_type = "patient" granularity
# Elixhauser score is derived from diagnoses_icd table

# compare new and old Elixhauser scores
data <- load_concepts(c("death", "bmi_all", "elix_wide", "age"), "miiv", 
                      id_type = "patient")
data <- replace_na(data, val = FALSE, vars = "death")
data[, c(index_var(data)) := NULL]
data[, Obesity := NULL]
data <- data[complete.cases(data)]

# fit a simple logistic model to understand the effect of comorbidities & bmi
logreg <- glm(death ~ . - subject_id - elix_wide - bmi_all + (bmi_all > 25),
              family = "binomial", data = data)
summary(logreg)

logreg_red <- glm(death ~ (bmi_all > 25) + elix_wide,
                  family = "binomial", data = data)
summary(logreg_red)

# frequency of comorbidities
cmbs <- names(data)[-(1:5)]
freq <- data.frame(condition = cmbs, cnt = colMeans(data[, cmbs, with = FALSE]),
                   coef = coef(logreg)[-c(1, 2, 32)])

cowplot::plot_grid(
  ggplot(freq, aes(x = factor(condition), y = cnt)) +
    theme_minimal() + geom_col() +
    # Other geoms and aesthetics can be added here
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)),
  ggplot(freq, aes(x = factor(condition), y = coef)) +
    theme_minimal() + geom_col() +
    # Other geoms and aesthetics can be added here
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)),
  ncol = 1L
)

# hypertension / diabetes + protective effect

ggplot(
  data, aes(x = bmi_all, fill = factor(DM))
) + geom_density(alpha = 0.3) + theme_minimal() +
  xlim(c(0, 50))

# do these negative coefficients disappear once you condition age?
# => yes, they do (partly)

# tipping over hypothesis is the most likely explanation for negative coefficients
# or association of HTN with obesity / fat mass / protective effect

# final task: get a top list of features
freq <- freq[order(freq$cnt, decreasing = TRUE), ]
freq <- freq[order(freq$coef, decreasing = TRUE), ]


# make a list of 5-10 features
# criterion: all with coef > 0.2, all with prevalence > 0.1 and coef > 0.05
cmb_sel <- freq[freq$coef > 0.2 | (freq$cnt > 0.1 & freq$coef > 0.1), ]$condition
# how many features is this?

# compute fi_crit
fi_crit <- function(r) {
  
  cor_mat <- t(r) %*% r / nrow(r)
  fi_crit <- min(
    min(1 - diag(cor_mat)),
    min(1 - sqrt(cor_mat[upper.tri(cor_mat)]))
  )
  
  fi_crit
}

fi_crit(as.matrix(data[, cmb_sel, with=FALSE]))


# check the coefficients / prevalence on ANZICS
anz_cmbs <- c("chr_resp", "chr_cvs", "chr_liv", "chr_ren", "chr_imun", 
              "chr_imunrx", "aids", "hep_fail", "lymphoma", "met_canc", 
              "leukaemia", "chr_imunsup", "cirrhosis")

anz_dat <- load_concepts(c(anz_cmbs, "site", "bmi", "age"), "anzics")
anz_dat <- replace_na(anz_dat, 65, vars = "age")
anz_dat <- anz_data[!is.na(bmi)]
anz_dat <- replace_na(anz_dat, 0)
anz_dat[,  anz_score := rowSums(.SD), .SDcols = anz_cmbs]

anz_freq <- anz_dat[, mean(anz_score), by = "site"]
anz_freq$site <- factor(anz_freq$site, 
                        levels = anz_freq$site[order(anz_freq$V1, decreasing = TRUE)])
ggplot(anz_freq, aes(x = site, y = V1)) +
  geom_col() + theme_minimal()

# first logistic model
anz_dat <- load_concepts(c(anz_cmbs, "site", "bmi", "age", "death"), "anzics")
anz_dat[, c(index_var(anz_dat)) := NULL]
anz_dat <- replace_na(anz_dat, FALSE, vars = "death")
anz_dat <- replace_na(anz_dat, 65, vars = "age")
anz_dat <- anz_dat[!is.na(bmi)]
anz_dat <- replace_na(anz_dat, 0)

# fit a logistic model and look at coefficients
anz_logreg <- glm(death ~ . - site - ICUStayID - bmi + (bmi > 25),
                  data = anz_dat, family = "binomial")

summary(anz_logreg)
