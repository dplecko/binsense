#' * three-step solution * 


#' * First step - load *
load_data <- function(src = "miiv", anzics_diag_subset = NULL, 
                      anzics_year_subset = NULL) {
  
  ind <- load_concepts(c("bmi", "death", "sex", "age"), src)
  ind[, c(index_var(ind)) := NULL]
  ind <- ind[!is.na(bmi) & !is.na(sex)]
  ind[is.na(death), death := FALSE]
  
  if (src != "anzics") {
    
    ils <- load_concepts("sofa", src, explicit_wins = hours(24L), 
                         keep_components = TRUE)
    ils[, c(index_var(ils), "sofa") := NULL]
    ind <- merge(ind, ils, all.x = TRUE)
    ind <- replace_na(ind, val = 0, type = "const")
  } else {
    ils <- load_concepts("apache_iii_risk", src)
    ind <- merge(ind, ils, all.x = TRUE)
    risk_med <- median(ind$apache_iii_risk, na.rm = TRUE) # compute median
    ind[is.na(apache_iii_risk), apache_iii_risk := risk_med] # impute median
    
    if (!is.null(anzics_diag_subset)) {
      diag <- load_concepts("apache_iii_diag", src)
      ind <- merge(ind, diag)
      ind <- ind[apache_iii_diag %in% anzics_diag_subset]
      ind[, apache_iii_diag := NULL]
    }
    
    if (!is.null(anzics_year_subset)) {
      year <- load_concepts("adm_year", src)
      ind <- merge(ind, year)
      ind <- ind[adm_year %in% anzics_year_subset]
      ind[, adm_year := NULL]
    }
  }
  
  if (src %in% c("mimic_demo", "mimic", "miiv", "aumc")) {
    dgs <- load_concepts("diag", src)
    ind <- merge(ind, dgs, all.x = TRUE)
  }
  
  
  if (any(is.na(ind))) browser()
  
  ind[!complete.cases(ind)]
  
  breaks <- c(0, 18, 25, 30, 35, Inf)
  # tags <- paste("<", breaks[-1])
  ind[, bmi_bin := factor(.bincode(bmi, breaks),
                          labels = c("<18", "<25", "<30", "<35", ">35"))]

  ind$bmi_bin <- relevel(ind$bmi_bin, ref = "<25")
  ind$sex <- factor(ind$sex)
  ind$sex <- relevel(ind$sex, ref = "Male")
  
  ind[, c("bmi", id_vars(ind)) := NULL]
  ind
}


#' * Second step - fit logistic models *
fit_log <- function(dat, type = c("univar", "multivar")) {
  
  type <- match.arg(type, c("univar", "multivar"))
  
  if (type == "univar") {
    mod <- glm(death ~ bmi_bin * sex, data = dat)
  } else {
    mod <- glm(death ~ bmi_bin * sex + ., data = dat)
  }
  
  mod
}


#' * Third step - plot the logistic coefficients * 
plot_or <- function(mod) {
  
  df <- summary(mod)$coefficients
  df <- rbind(df, `bmi_bin<25`=0)
  df <- rbind(df, `bmi_bin<25:sexFemale`=0)
  df <- data.table(df, keep.rownames = TRUE)
  
  dif_impact <- df[rn == "sexFemale"]$Estimate
  
  df <- df[grepl("bmi_bin", rn)]
  df[, group := ifelse(grepl("Female", rn), "Female", "Male")]
  
  df <- setorderv(df, c("group", "rn"))
  df[group == "Female"]$Estimate <- df[group == "Female"]$Estimate + 
    df[group == "Male"]$Estimate + dif_impact
  
  df[group == "Female", Estimate := Estimate + dif_impact]
  df[, bmi_bin := gsub("bmi_bin|:.*", "", df$rn)]
  df$bmi_bin <- factor(df$bmi_bin)
  
  ggplot(df, aes(x = as.integer(bmi_bin), y = exp(Estimate), color = group)) +
    geom_point() + geom_line() + theme_bw() +
    geom_errorbar(aes(ymin = exp(Estimate - 1.96 * `Std. Error`),
                      ymax = exp(Estimate + 1.96 * `Std. Error`)),
                  linewidth = 0.75, width = 0.25) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    scale_x_continuous(labels = levels(df$bmi_bin), 
                       breaks = seq_along(levels(df$bmi_bin))) +
    ylab("Estimate Odds Ratio") + xlab("BMI group") + 
    theme(
      legend.position = c(0.8, 0.8),
      legend.box.background = element_rect()
    ) + scale_color_discrete(name = "Sex")

}
