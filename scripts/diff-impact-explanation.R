
#' * explaining variations approach *
chr_il <- c("chr_resp", "chr_cvs", "chr_liv", "chr_ren", "chr_imun", "chr_imunrx", 
            "aids", "hep_fail", "lymphoma", "met_canc", "leukaemia", "chr_imunsup", 
            "cirrhosis")

z_set <- c(chr_il, "bmi_bin", "age")

w_set <- c("apache_iii_diag", "apache_iii_risk")

dat <- load_concepts(c("death", "sex", z_set, w_set), "anzics")
dat <- dat[complete.cases(dat)]
dat[, apache_iii_diag := as.factor(apache_iii_diag)]
dat$death <- as.integer(dat$death)

anz_expl <- faircause::fairness_cookbook(
  data = dat, X = "sex", W = w_set, Z = z_set, Y = "death", x0 = "Female", 
  x1 = "Male"
)

autoplot(anz_expl, decompose = "xspec")

#' * logistic regression exploration *

mod0 <- glm(death ~ sex, dat,
           family = "binomial")
summary(mod0)

mod <- glm(death ~ apache_iii_risk + factor(apache_iii_diag) + sex, dat,
           family = "binomial")
summary(mod)

ggsave(file.path(root, "results", "tv-decomposed.png"), width = 6, height = 5)

library(ricu)

dat <- load_concepts(c("death", "sex", "bmi"), "anzics")

dat[, mean(death), by = "sex"]

dat[is.na(death), death := FALSE]

for (col in chr_il) {
  dat[is.na(get(col)), c(col) := 0]
}


chron <- dat[, lapply(.SD, mean), by = "sex", .SDcols = chr_il]


# not chronic illness

plot(as.numeric(chron[1, -1]), lty = 2, pch = 19)
lines(as.numeric(chron[2, -1]), lty = 2, pch = 19, col = "red")
# ggplot(
#   melt([c(1, 2)]), 
#   aes(x = variable, y = value, fill = sex)) + 
#   geom_line(position = "identity", alpha = 0.4)


ggplot(dat, aes(x = factor(apache_iii_diag), fill = sex, y = death)) +  
  geom_bar(aes(y = after_stat(count)/sum(after_stat(count)))) + 
  scale_y_continuous(labels = scales::percent)

diag <-
  dat[!is.na(sex), list(mean(death), mean(sex == "Male"), 122 * .N / nrow(dat[!is.na(sex)])), by = "apache_iii_diag"]

diag <- setorderv(diag, cols = c("V1"), order = -1L)
diag[, diag_ord := seq_len(nrow(diag))]
diag[, V3 := V3 / 10]

ggplot(
  melt(diag, id.vars = c("apache_iii_diag", "diag_ord")), 
  aes(x = diag_ord, y = value, color = variable)) + 
  geom_line() + theme_bw() +
  geom_hline(yintercept = 0.57, linetype = "dashed", color = "red", linewidth = 2) +
  scale_y_continuous(labels = scales::percent)

#' * APACHE-III risk density *
dat_age <- load_concepts(c("age", "sex", "apache_iii_risk"), "anzics")
dat_age <- dat_age[age > 55]
dat_age <- dat_age[complete.cases(dat_age)]
ggplot(dat_age, aes(x = apache_iii_risk, fill = sex, color = sex)) +
  #geom_density(alpha = 0.4) + 
  stat_ecdf(geom = "step") +
  theme_bw() 

ggsave(file.path(root, "results", "AP3_rod.png"), width = 5, height = 3)

dat[, mean(apache_iii_risk), by = "sex"]

#' * Age ecdf *
dat_age <- load_concepts(c("age", "sex", "apache_iii_diag"), "anzics")
dat_age <- dat_age[!(apache_iii_diag %in% c(1800+1:3, 902, 903, 1408, 2201, 704, 1903, 1904))]
dat_age <- dat_age[complete.cases(dat_age)]
ggplot(dat_age, aes(x = age, fill = sex, color = sex)) +
  #geom_density(alpha = 0.4) + 
  stat_ecdf(geom = "step") +
  theme_bw()

ggsave(file.path(root, "results", "age_miiv.png"), width = 5, height = 3)


#' * Age-driving diagnoses *
dat <- load_concepts(c("age", "sex", "apache_iii_diag","death"), "anzics")
dat[is.na(death), death := FALSE]
dat <- dat[!(apache_iii_diag %in% c(1800+1:3, 902, 903, 1408, 2201, 704, 1903, 1904))]
dat <- dat[complete.cases(dat)]

dat <- dat[age <= 50 & age >= 25, .N, by = c("apache_iii_diag", "sex")]

dat <- split(dat[, -2], dat$sex)
dat <- merge(dat[[1]], dat[[2]], by = "apache_iii_diag", all = TRUE)
dat[, gain := N.x - N.y]

setorderv(dat, "gain", 1L)
tail(dat, n = 20L)

plot(dat[-1]$N.x, pch = 19, col = "red")
points(dat[-1]$N.y, col = "blue")


