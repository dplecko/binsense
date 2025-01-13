
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

elixhauser_cb <- function(sex, ...) {
  
  icd_to_elix <- function(x, ...) {
    
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
    w_all
  }
  
  id_type <- id_var(sex)
  diag <- load_src("diagnoses_icd", "miiv")
  id_col <- which(names(diag) == id_type)
  x <- as_id_tbl(diag, id_vars = id_col)
  res <- miiv_elix_wide2(x, ...)
  res[, elixhauser := sum(.SD), 
      .SDcols = setdiff(names(res), id_vars(res))]
  res
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

wide_elix <- function(elix_all, ...) {
  
  res <- merge(elix_all, get_elix())
  res <- rename_cols(res, "elix_wide", "elix_all")
  res
}

get_elix <- function(src = "miiv") {
  
  load(file.path(root, "data", "miiv_elix_wide.rda"))
  w_all
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