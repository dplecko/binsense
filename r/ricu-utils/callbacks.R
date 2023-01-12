
ismv_cb <- function(mech_vent, ...) {

  mech_vent <- mech_vent[, .N, by = c(id_vars(mech_vent))]
  mech_vent[, is_mv := N > 0]
  mech_vent[, N := NULL]
}

anzics_binary <- function(x, val_var, ...) {
  
  x[, c(val_var) := ifelse(get(val_var) == 1, 1, 0)]
  x
}

anzics_death <- function(x, ...) {
  
  x[, DIED := ifelse(DIED == 1, TRUE, FALSE)]
  
  x
}

anzics_sex <- function(x, ...) {

  x[, SEX := ifelse(SEX == "M", "Male", "Female")]
  
  x
}

eicu_diag <- function(...) {
  browser()
}

bmi_binary <- function(bmi, ...) {
  
  breaks <- c(-Inf, 25, Inf)
  
  bmi[, bmi_bin := factor(.bincode(bmi, breaks))]
  levels(bmi[["bmi_bin"]]) <- c("<25", ">25")
  
  id_var <- id_vars(bmi)
  bmi[, c(id_var, "bmi_bin"), with = F]
  
}

charlson_callback <- function(x, ...) {

  sub_var <- setdiff(names(x), meta_vars(x))
  if (sub_var == "icd9code") {
    
    x[, c(sub_var) := gsub(",.*", "", get(sub_var))]
    
  }
  
  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  intm <- rowSums(comorbid_charlson(intm))
  
  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )
  
  names(res) <- names(x)
  
  res
}

elix_callback <- function(x, ...) {
  
  sub_var <- setdiff(names(x), meta_vars(x))
  if (sub_var == "icd9code") {
    
    x[, c(sub_var) := gsub(",.*", "", get(sub_var))]
    
  }
  
  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  
  rm_cols <- c("DM", "DMcx", "Obesity", "WeightLoss")
  intm <- comorbid_elix(intm)
  rm_cols <- which(colnames(intm) %in% rm_cols)
  intm <- rowSums(intm[, -rm_cols])
  #browser()
  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )
  
  names(res) <- names(x)
  
  res
}

elix_miiv <- function(x, group_var, ...) {
  
  sub_var <- setdiff(names(x), c(meta_vars(x), group_var))
  
  get_elix <- function(x, sub_var) {
    intm <- data.frame(
      pid = id_col(x),
      icd = x[[sub_var]]
    )
    
    rm_cols <- c("DM", "DMcx", "Obesity", "WeightLoss")
    intm <- comorbid_elix(intm)
    rm_cols <- which(colnames(intm) %in% rm_cols)
    intm <- rowSums(intm[, -rm_cols])
    res <- id_tbl(
      id = as.integer(names(intm)),
      val = intm, id_vars = "id"
    )
    
    names(res) <- c(id_var(x), "icd_code")
    
    res
  }
  # browser()
  rbind(
    get_elix(x[get(group_var) == 9], sub_var),
    get_elix(x[get(group_var) == 10], sub_var)
  )
}

icd10_callback <- function(x, group_var, val_var, ...) {

  x[icd_version == 10, c(id_vars(x), "icd_code"), with=FALSE]
  
}

icd9_callback <- function(x, group_var, val_var, ...) {
  x1 <- x[icd_version == 9, c(id_vars(x), "icd_code"), with=FALSE]
  x1[, version := 9]
  
  x2 <- x[icd_version == 10, c(id_vars(x), "icd_code"), with=FALSE]
  x2[, version := 10]
  icd <- rbind(x1, x2)
  save(icd, file = file.path(root, "data", "icd.rda"))
  browser()
  res
}

icd_merge_callback <- function(icd9, icd10, ...) {
  browser()
  icd9[, version := 9]
  icd9 <- rename_cols(icd9, "icd_code", "icd9")
  icd10[, version := 10]
  icd10 <- rename_cols(icd10, "icd_code", "icd10")
  rbind(icd9, icd10)
  
}

elix_miiv_full <- function(x, group_var, ...) {
  
  sub_var <- setdiff(names(x), c(meta_vars(x), group_var))
  
  get_elix_full <- function(x, sub_var) {
    intm <- data.frame(
      pid = id_col(x),
      icd = x[[sub_var]]
    )
    
    rm_cols <- c("DM", "DMcx", "Obesity", "WeightLoss")
    intm <- comorbid_elix(intm)
    rm_cols <- which(colnames(intm) %in% rm_cols)
    intm <- intm[, -rm_cols]
    # res <- id_tbl(
    #   id = as.integer(names(intm)),
    #   val = intm, id_vars = "id"
    # )
    # 
    # names(res) <- c(id_var(x), "icd_code")
    
    intm
  }
  browser()
  rbind(
    get_elix_full(x[get(group_var) == 9], sub_var),
    get_elix_full(x[get(group_var) == 10], sub_var)
  )
}

DM910_callback <- function(x, val_var, ...) {
  
  if (val_var == "icd9code") {
    
    x[, c(val_var) := gsub(",.*", "", get(val_var))]
    
  }
  
  DM_map <- list(
    DM = c(icd9_map_charlson$DM, icd10_map_charlson$DM),
    DMcx = c(icd9_map_charlson$DMcx, icd10_map_charlson$DMcx)
  )
  
  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  intm <- rowSums(comorbid(intm, map = DM_map)[, c("DM", "DMcx")]) > 0
  
  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )
  
  names(res) <- names(x)
  
  res
}

mimic_adm_diag <- function(x, val_var, ...) {
  
  mapp <- list(
    OTH = c("DENT", "PSYCH", "NBB", "NB", "ENT", "GU", "PSURG"),
    GYN = c("GYN", "OBS")
  )
  for (i in seq_along(mapp)) {
    x[get(val_var) %in% mapp[[i]], c(val_var) := names(mapp)[i]]
  }
  x
}

tw_avg_gluc <- function(glu, icu_end, upto = hours(24L), 
                        hypo_censoring = TRUE, ...) {
  
  ind <- index_var(glu)
  limits <- merge(glu[, list(first_obs = min(get(ind))), 
                      by = c(id_vars(glu))], 
                  icu_end)
  limits[, start := pmin(hours(0L), first_obs)]
  if (hypo_censoring) {
    hg <- glu[, list(hypo_time = head(get(ind)[glu <= 70], 1L)), 
              by = c(id_vars(glu))]
    limits <- merge(limits, hg, all.x = TRUE)
    limits[is.na(hypo_time), hypo_time := hours(Inf)]
    limits[, end := pmin(icu_end, hypo_time - hours(1L))]
  }
  
  glu <- fill_gaps(glu, limits = limits)
  glu[, glu := data.table::nafill(glu, "locf"), by = eval(id_vars(glu))]
  
  glu[get(ind) >= 0L & get(ind) <= upto]
  
  glu[, list(tw_avg_glu = mean(glu, na.rm = T)), by = eval(id_vars(glu))]
  
}

is_hypo_cb <- function(glu, interval, ...) {
  idx <- index_var(glu)
  glu[get(idx) <= hours(24L), list(is_hypo = any(glu <= 70)), 
      by = c(id_vars(glu))]
  
}

SMK_callback <- function(x, val_var, ...) {
  
  if (val_var == "icd9code") {
    
    x[, c(val_var) := gsub(",.*", "", get(val_var))]
    
  }
  
  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  intm <- rowSums(
    icd9_comorbid(
      intm, map = list(
        Smoking = c("3051", "305.1"))
      )[, c("Smoking"), drop = FALSE]) > 0
  
  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )
  
  names(res) <- names(x)
  
  res
}

SMK_miiv <- function(x, group_var, ...) {
  
  sub_var <- setdiff(names(x), c(meta_vars(x), group_var))
  
  get_elix <- function(x, sub_var, icd_fn) {
    intm <- data.frame(
      pid = id_col(x),
      icd = x[[sub_var]]
    )
    
    intm <- rowSums(
      icd_fn(intm, map = list(
        Smoking = c("3051", "305.1", "Z72.0", "Z87.891", "F17.2"))
      )[, c("Smoking"), drop = FALSE]) > 0
    res <- id_tbl(
      id = as.integer(names(intm)),
      val = intm, id_vars = "id"
    )
    
    names(res) <- c(id_var(x), "icd_code")
    
    res
  }
  
  rbind(
    get_elix(x[get(group_var) == 9], sub_var, icd9_comorbid),
    get_elix(x[get(group_var) == 10], sub_var, icd10_comorbid)
  )
}

get_icd <- function(src = "miiv", index = "elix") {
  
  load(file.path(root, "data", "icd.rda"))
  
  get_elix <- function(x, sub_var, index) {
    intm <- data.frame(
      pid = id_col(x),
      icd = x[[sub_var]]
    )
    
    rm_cols <- c("Obesity", "WeightLoss")
    if (index == "charlson") {
      intm <- comorbid_charlson(intm)
    } else {
      intm <- comorbid_elix(intm)
    }
    rm_cols <- which(colnames(intm) %in% rm_cols)
    if (length(rm_cols) > 0) intm <- intm[, -rm_cols]
    res <- id_tbl(
      id = as.integer(rownames(intm)), id_vars = "id"
    )
    
    cbind(res, intm)
  }
  
  get_elix(icd, "icd_code", index)
  
}

expit <- function(x) {
  exp(x) / (1 + exp(x))
}
