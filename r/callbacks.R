
weight_anzics_cb <- function(x, ...) {
  
  height <- x$HEIGHT
  weight <- x$WEIGHT
  age <- x$AGE
  
  height2 <- height
  weight2 <- weight
  weight[which(height2 < 100)] <- height2[which(height2 < 100)]
  height[which(height2 < 100)] <- weight2[which(height2 < 100)]
  weight2[which(weight2 == 0)] <- NA
  height2[which(weight2 == 0)] <- NA
  height[which(height == 0)] <- NA
  height[which(height > 250)] <- NA
  height[which(height < 30)] <- NA
  height[which(height < 50 & age > 5)] <- NA
  height[which(height < 100 & age > 10)] <- NA
  weight[which(weight == 0)] <- NA
  weight[which(weight > 350)] <- NA
  weight[which(weight > 200 & height < 120)] <- NA
  weight[which(weight < 20 & age > 10)] <- NA
  weight[which(weight < 10 & age > 5)] <- NA
  
  x[, WEIGHT := weight]
  
  x[, c(id_vars(x), "WEIGHT"), with=FALSE]
}

height_anzics_cb <- function(x, ...) {
  
  height <- x$HEIGHT
  weight <- x$WEIGHT
  age <- x$AGE
  
  height2 <- height
  weight2 <- weight
  weight[which(height2 < 100)] <- height2[which(height2 < 100)]
  height[which(height2 < 100)] <- weight2[which(height2 < 100)]
  weight2[which(weight2 == 0)] <- NA
  height2[which(weight2 == 0)] <- NA
  height[which(height == 0)] <- NA
  height[which(height > 250)] <- NA
  height[which(height < 30)] <- NA
  height[which(height < 50 & age > 5)] <- NA
  height[which(height < 100 & age > 10)] <- NA
  weight[which(weight == 0)] <- NA
  weight[which(weight > 350)] <- NA
  weight[which(weight > 200 & height < 120)] <- NA
  weight[which(weight < 20 & age > 10)] <- NA
  weight[which(weight < 10 & age > 5)] <- NA
  
  x[, HEIGHT := height]
  
  x[, c(id_vars(x), "HEIGHT"), with=FALSE]
}

anzics_adm_type_cb <- function(adm, elective, ...) {
  
  adm_type <- function(adm, elective) {
    ifelse(
      adm == "med", "med", ifelse(elective, "elective_surgery", "emergency_surgery")
    )
  }
  res <- merge(adm, elective, all = TRUE)
  res[, adm_type := adm_type(adm, elective)]
  res[, c(id_vars(res), "adm_type"), with=FALSE]
}

is_chr_cb <- function(cmbd_anzics, ...) {
  
  cmbd_anzics[, is_chr := cmbd_anzics > 0]
  cmbd_anzics[, c(id_vars(cmbd_anzics), "is_chr"), with=FALSE]
}

is_vent_cb <- function(...) {
  
  res <- Reduce(function(x, y) merge(x, y, all = TRUE), list(...)[1:3])
  res[, is_vent := is_invasive > 0 | is_invasive2 > 0 | is_noninvasive > 0]
  res[is.na(is_vent), is_vent := FALSE]
  res[, c(id_vars(res), "is_vent"), with=FALSE]
}

anzics_diab_cb <- function(x, ...) {
  
  x[, DIABETES := DIABETES != 5]
  x
}

cmbd_anzics_cb <- function(...) {
  
  res <- Reduce(function(x, y) merge(x, y, all = TRUE), list(...)[1:13])
  res <- replace_na(res, 0)
  res[, cmbd_anzics := rowSums(res[, -c(id_vars(res)), with=FALSE])]
  res[, c(id_vars(res), "cmbd_anzics"), with=FALSE]
}

miiv_elix_dir <- function(x, ...) {
  
  ch9 <- icd9_comorbid_elix(x[icd_version == 9])
  ch10 <- icd10_comorbid_elix(x[icd_version == 10])
  
  make_long <- function(x, id_name) {
    
    res <- id_tbl(id = as.integer(rownames(x)))
    res <- cbind(res, x)
    res <- rename_cols(res, id_name, "id")
    res <- melt.data.table(res, id.vars = id_name)
    as_id_tbl(res)
  }

  ch <- rbind(
    make_long(ch9, id_vars(x)),
    make_long(ch10, id_vars(x))
  )
  browser()
  ch <- ch[!(variable %in% c("Obesity", "WeightLoss"))]
  ch <- ch[, list(cmb = max(value, na.rm = TRUE)), by = c(id_vars(ch), "variable")]
  ch[, list(icd_code = sum(cmb, na.rm = TRUE)), by = c(id_vars(ch))]
}

miiv_elix_wide <- function(x, ...) {
  
  ch9 <- icd9_comorbid_elix(x[icd_version == 9])
  ch10 <- icd10_comorbid_elix(x[icd_version == 10])
  
  make_wide <- function(x, id_name) {
    
    res <- id_tbl(id = as.integer(rownames(x)))
    res <- cbind(res, x)
    res <- rename_cols(res, id_name, "id")
    res
  }
  
  w9 <- make_wide(ch9, id_vars(x))
  w10 <- make_wide(ch10, id_vars(x))
  
  w_all <- rbind(w9, w10)
  w_all <- w_all[, lapply(.SD, max), by = c(id_vars(w_all))]
  save(w_all, file = file.path(root, "data", "miiv_elix_wide.rda"))
  w_all[, icd_code := rowSums(w_all[, -c(id_vars(w_all)), with=FALSE])]
}

miiv_charlson_wide <- function(x, ...) {
  
  ch9 <- icd9_comorbid_charlson(x[icd_version == 9])
  ch10 <- icd10_comorbid_charlson(x[icd_version == 10])
  
  make_wide <- function(x, id_name) {
    
    res <- id_tbl(id = as.integer(rownames(x)))
    res <- cbind(res, x)
    res <- rename_cols(res, id_name, "id")
    res
  }
  
  w9 <- make_wide(ch9, id_vars(x))
  w10 <- make_wide(ch10, id_vars(x))
  
  w_all <- rbind(w9, w10)
  w_all <- w_all[, lapply(.SD, max), by = c(id_vars(w_all))]
  save(w_all, file = file.path(root, "data", "miiv_charlson_wide.rda"))
  w_all[, icd_code := rowSums(w_all[, -c(id_vars(w_all)), with=FALSE])]
}

wide_elix <- function(elix_all, ...) {
  
  res <- merge(elix_all, get_elix())
  res <- rename_cols(res, "elix_wide", "elix_all")
  res
}

get_elix <- function(src = "miiv") {
  
  load(file.path(root, "data", "miiv_elix_wide.rda"))
  w_all
}

get_charlson <- function(src = "miiv") {
  
  load(file.path(root, "data", "miiv_charlson_wide.rda"))
  w_all
}

miiv_charlson_dir <- function(x, ...) {

  ch9 <- icd9_comorbid_charlson(x[icd_version == 9])
  ch10 <- icd10_comorbid_charlson(x[icd_version == 10])
  
  make_long <- function(x, id_name) {
    
    res <- id_tbl(id = as.integer(rownames(x)))
    res <- cbind(res, x)
    res <- rename_cols(res, id_name, "id")
    res <- melt.data.table(res, id.vars = id_name)
    as_id_tbl(res)
  }
  
  ch <- rbind(
    make_long(ch9, id_vars(x)),
    make_long(ch10, id_vars(x))
  )
  
  ch <- ch[, list(cmb = max(value, na.rm = TRUE)), by = c(id_vars(ch), "variable")]
  ch[, list(icd_code = sum(cmb, na.rm = TRUE)), by = c(id_vars(ch))]
}

cmb_charlson_cb <- function(charlson9, charlson10, ...) {
  
  charlson9 <- rename_cols(charlson9, "charlson", "charlson9")
  charlson10 <- rename_cols(charlson10, "charlson", "charlson10")
  rbind(charlson9, charlson10)
}

miiv_charlson10_cb <- function(icd_codes, ...) {
  
  v10 <- icd10_comorbid_charlson(icd_codes)
  v10 <- id_tbl(
    subject_id = as.integer(rownames(v10)),
    icd_code = rowSums(v10)
  )
  
  v10
}

miiv_charlson9_cb <- function(icd_codes, ...) {

  v9 <- icd9_comorbid_charlson(icd_codes)
  v9 <- id_tbl(
    subject_id = as.integer(rownames(v9)),
    icd_code = rowSums(v9)
  )
  
  v9
}

miiv_adm_epi_cb <- function(x, ...) {
  
  x <- merge(x, list(...)$env$icustays[, c("stay_id", "intime")], by = "stay_id")
  x <- setorderv(x, cols = c("subject_id", "intime"))
  x[, adm_episode := seq_len(.N), by = "subject_id"]
  
  x <- x[, c(id_vars(x), "adm_episode"), with=FALSE]
  rename_cols(x, "subject_id", "adm_episode")
}

bmi_all_cb <- function(bmi, bmi_omr, ...) {
  
  bmi <- rename_cols(bmi, "bmi_all", "bmi")
  bmi_omr <- rename_cols(bmi_omr, "bmi_all", "bmi_omr")
  
  x <- rbind(bmi, bmi_omr)
  x[, list(bmi_all = mean(bmi_all)), by = c(id_vars(x))]
}

anzics_omr_bmi_cb <- function(x, ...) {
  
  x[integer(0L)]
}

miiv_omr_bmi_cb <- function(x, ...) {
  
  x[, result_value := as.numeric(result_value)]
  x[, bmi_omr := mean(result_value, na.rm = TRUE), by = c(id_vars(x))]
  
  x
}

los_hosp_anzics_cb <- function(x, ...) {
  
  x[, HOSP_HRS := HOSP_HRS / 24]
  x
}

los_icu_anzics_cb <- function(x, ...) {
  
  x[, ICU_HRS := ICU_HRS / 24]
  x
}

pci_or_death_callback <- function(pci, death, ...) {
  
  x <- merge(pci, death, all = TRUE)
  x[, pci_or_death := death | pci]
  x[is.na(pci_or_death), pci_or_death := FALSE]
  x[, c(id_vars(x), "pci_or_death"), with=FALSE]
}

pci_callback <- function(los_icu, ...) {
  
  los_icu[, pci := los_icu >= 10]
  los_icu[, c(id_vars(los_icu), "pci"), with=FALSE]
}

anzics_adm <- function(x, val_var, ...) {
  
  diag_to_adm <- function(x) ifelse(x < 1200, "med", "surg")
  
  x[, adm := diag_to_adm(get(val_var))]
  x[, c(val_var) := NULL]
  x <- setnames(x, "adm", val_var)
  x
}

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
  
  # browser()
  x$uid <- seq_len(nrow(x))
  x[icd_version == 10, c(id_vars(x), "icd_code", "uid"), with=FALSE]
  
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

bin_bmi <- function(bmi, ...) {
  
  breaks <- c(-Inf, config("bmi-bins")[["who"]], Inf)
  
  bmi[, bmi_bins := factor(.bincode(bmi, breaks))]
  levels(bmi[["bmi_bins"]]) <- bin_labels(config("bmi-bins")[["who"]], "$kg/m^2$")
  
  id_var <- id_vars(bmi)
  bmi[, c(id_var, "bmi_bins"), with = FALSE]
  
}
