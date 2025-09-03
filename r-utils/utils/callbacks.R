
is_vent_cb <- function(...) {
  
  res <- Reduce(function(x, y) merge(x, y, all = TRUE), list(...)[1:3])
  res[, is_vent := is_invasive > 0 | is_invasive2 > 0 | is_noninvasive > 0]
  res[is.na(is_vent), is_vent := FALSE]
  res[, c(id_vars(res), "is_vent"), with=FALSE]
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
  
  browser()
  id_type <- id_var(sex)
  diag <- load_src("diagnoses_icd", "miiv")
  id_col <- which(names(diag) == id_type)
  x <- as_id_tbl(diag, id_vars = id_col)
  res <- miiv_elix_wide2(x, ...)
  res[, elixhauser := sum(.SD), 
      .SDcols = setdiff(names(res), id_vars(res))]
  res
}

# miiv_elix_wide <- function(x, ...) {
#   
#   browser()
#   
#   ch9 <- icd9_comorbid_elix(x[icd_version == 9])
#   ch10 <- icd10_comorbid_elix(x[icd_version == 10])
#   
#   make_wide <- function(x, id_name) {
#     
#     res <- id_tbl(id = as.integer(rownames(x)))
#     res <- cbind(res, x)
#     res <- rename_cols(res, id_name, "id")
#     res
#   }
#   
#   w9 <- make_wide(ch9, id_vars(x))
#   w10 <- make_wide(ch10, id_vars(x))
#   
#   w_all <- rbind(w9, w10)
#   w_all <- w_all[, lapply(.SD, max), by = c(id_vars(w_all))]
#   save(w_all, file = file.path(root, "data", "miiv_elix_wide.rda"))
#   w_all[, icd_code := rowSums(w_all[, -c(id_vars(w_all)), with=FALSE])]
# }

# 
miiv_elix_wide <- function(x, ...) {

  ch9 <- icd9_comorbid_elix(x[icd_version == 9])
  ch10 <- icd10_comorbid_elix(x[icd_version == 10])

  w_all <- rbind(ch9, ch10)
  w_all <- w_all[, lapply(.SD, max), by = c(id_vars(w_all))]
  save(w_all, file = file.path(root, "data", "miiv2_elix_wide.rda"))
  w_all[, icd_code := rowSums(w_all[, -c(id_vars(w_all)), with=FALSE])]
}

wide_elix <- function(elix_all, ...) {
  
  res <- merge(elix_all, get_elix())
  res <- rename_cols(res, "elix_wide", "elix_all")
  res
}

get_elix <- function(src = "miiv") {
  
  load(file.path(root, "data", "miiv2_elix_wide.rda"))
  w_all
}

bmi_all_cb <- function(bmi, bmi_omr, ...) {
  
  bmi <- rename_cols(bmi, "bmi_all", "bmi")
  bmi_omr <- rename_cols(bmi_omr, "bmi_all", "bmi_omr")
  
  x <- rbind(bmi, bmi_omr)
  x[, list(bmi_all = mean(bmi_all)), by = c(id_vars(x))]
}

miiv_omr_bmi_cb <- function(x, ...) {
  
  x[, result_value := as.numeric(result_value)]
  x[, bmi_omr := mean(result_value, na.rm = TRUE), by = c(id_vars(x))]
  
  x
}

ismv_cb <- function(mech_vent, ...) {

  mech_vent <- mech_vent[, .N, by = c(id_vars(mech_vent))]
  mech_vent[, is_mv := N > 0]
  mech_vent[, N := NULL]
}

bmi_binary <- function(bmi, ...) {
  
  breaks <- c(-Inf, 25, Inf)
  
  bmi[, bmi_bin := factor(.bincode(bmi, breaks))]
  levels(bmi[["bmi_bin"]]) <- c("<25", ">25")
  
  id_var <- id_vars(bmi)
  bmi[, c(id_var, "bmi_bin"), with = F]
  
}

# manual Elixhauser comorbidities

# Normalize function
# Normalize ICD code: remove dot, uppercase
norm_code <- function(dt) {
  dt[, icd_code := toupper(gsub("\\.", "", icd_code))]
}

# Compute comorbidity flags by stay_id
compute_flags <- function(dt, map) {
  out_list <- vector("list", length(map))
  names(out_list) <- names(map)
  
  for (cat in names(map)) {
    codes <- map[[cat]]
    tmp <- dt[icd_code %in% codes, .(flag = 1L), by = c(id_vars(dt))]
    if (nrow(tmp) == 0L) {
      out_list[[cat]] <- unique(dt[, c(id_vars(dt)), with=FALSE])[, (cat) := 0L]
    } else {
      setnames(tmp, "flag", cat)
      out_list[[cat]] <- tmp
    }
  }
  
  res <- Reduce(function(a, b) merge(a, b, by = id_vars(dt), all = TRUE), out_list)
  for (j in setdiff(names(res), id_vars(dt))) {
    set(res, which(is.na(res[[j]])), j, 0L)
  }
  
  res
}

# ICD-9 wrapper
icd9_comorbid_elix <- function(dt) {
  dt <- copy(dt)
  norm_code(dt)
  compute_flags(dt, icd9_elix_map)
}

# ICD-10 wrapper
icd10_comorbid_elix <- function(dt) {
  dt <- copy(dt)
  norm_code(dt)
  compute_flags(dt, icd10_elix_map)
}

icd9_elix_map <- list(
  chf = c("39891", "40211", "40291", "40411", "40413", "40491", "40493", paste0("428", 0:9)),
  
  arrhythmia = c(
    "42610", "42611", "42613", "4262", "4263", "4264", "42650", "42651", "42652", "42653",
    "42660", "42681", "42689", "4270", "4271", "4272", "42731", "42760", "4279", "7850", 
    "V450", "V533"
  ),
  
  valve = c(
    "3940", "3941", "3942", "3949", "3950", "3951", "3952", "3959", "3960", "3961", "3962", 
    "3963", "3968", "3969", "3970", "3971", "4240", "4241", "4242", "4243", "42490", "42491",
    "7463", "7464", "7465", "7466", "V422", "V433"
  ),
  
  pulm.circ = c(paste0("416", 0:9), "4179"),
  
  pvd = c(
    paste0("440", 0:9), "4412", "4414", "4417", "4419", paste0("443", 1:9), 
    "4471", "5571", "5579", "V434"
  ),
  
  htn = c("4011", "4019"),
  
  htncx = c("40210", "40290", "40410", "40490", "40511", "40519", "40591", "40599"),
  
  paralysis = c(
    paste0("342", 0:9), "3430", "3431", "3432", "3433", "3434", "3438", "3439",
    "3440", "3441", "3442", "3443", "3444", "3445", "3446", "3449"
  ),
  
  neuro.other = c(
    "3319", "3320", "3334", "3335", paste0("334", 0:9), paste0("335", 0:9), "340", 
    "3411", "3418", "3419", sprintf("345%02d", 0:9), "3453", "3454", "3455", "3459", 
    "3481", "3483", "7803", "7843"
  ),
  
  chronic.pulm = c(
    "490", "4910", "4911", "4912", "4918", "4919", "4920", "4928", 
    "49300", "49310", "49320", "49390", "494", paste0("495", 0:9), 
    "496", "500", "501", "502", "503", "504", "505", "5064"
  ),
  
  dm.uncomp = sprintf("250%02d", 0:33),
  dm.comp = sprintf("250%02d", c(40:73, 90:93)),
  
  hypothyroid = c("243", "2440", "2441", "2442", "2448", "2449"),
  
  renal = c("40311", "40391", "40412", "40492", "585", "586", "V420", "V451", "V560", "V568"),
  
  liver = c(
    "07032", "07033", "07054", "4560", "4561", "45620", "45621", "5710", "5712", "5713", 
    "57140", "57141", "57142", "57149", "5715", "5716", "5718", "5719", "5723", "5728", 
    "V427"
  ),
  
  pud = c("53170", "53190", "53270", "53290", "53370", "53390", "53470", "53490", "V1271"),
  
  hiv = c("042", "043", "044"),
  
  lymphoma = c(
    sprintf("200%02d", 0:9), "20190", "20200", "20210", "20220", "20280", "20300", 
    "20301", "20380", "20381", "2386", "2733", "V1071", "V1072", "V1079"
  ),
  
  mets = c("1960", "1961", "1962", "1963", "1965", "1966", "1968", "1969", 
           "1970", "1971", "1972", "1973", "1974", "1975", "1976", "1977", 
           "1980", "1981", "1982", "1983", "1984", "1985", "1986", "1987", 
           "19881", "19882", "19889", "1990", "1991"),
  
  solid.tumor = c(
    sprintf("%03d", 140:172), sprintf("174%01d", 0:9), sprintf("175%01d", 0:9),
    sprintf("179", 0), sprintf("180%01d", 0:9), sprintf("181", 0),
    sprintf("182%01d", 0:9), sprintf("183%01d", 0:9), sprintf("184%01d", 0:9),
    sprintf("185", 0), sprintf("186%01d", 0:9), sprintf("187%01d", 0:9),
    sprintf("188%01d", 0:9), sprintf("189%01d", 0:9)
  ),
  
  rheum = c("7010", paste0("710", 0:9), paste0("714", 0:9), paste0("720", 0:9), "725"),
  
  coag = c(paste0("286", 0:9), "2871", "2873", "2874", "2875"),
  
  obesity = "2780",
  
  wt.loss = c("260", "261", "262", "2630", "2631", "2632", "2638", "2639"),
  
  lytes = paste0("276", 0:9),
  
  anemia.loss = "2800",
  
  anemia.def = c(paste0("280", 1:9), paste0("281", 0:9), "2859"),
  
  etoh = c("2911", "2912", "2915", "2918", "2919", paste0("3039", 0:3), paste0("3050", 0:3), "V113"),
  
  drugs = c(
    "2920", paste0("2928", 2:9), "2929", paste0("304", sprintf("%02d", 0:93)), 
    paste0("3052", 0:3), paste0("3053", 0:3), paste0("3054", 0:3), paste0("3055", 0:3), 
    paste0("3056", 0:3), paste0("3057", 0:3), paste0("3059", 0:3)
  ),
  
  psychoses = c(paste0("295", sprintf("%02d", 0:9)), paste0("296", sprintf("%02d", 0:9)),
                paste0("297", sprintf("%02d", 0:9)), paste0("298", sprintf("%02d", 0:9))),
  
  depression = c("3004", "30112", "3090", "3091", "311")
)


icd10_elix_map <- list(
  chf = c("I099", "I110", "I130", "I132", "I255", "I420", "I425", "I426", "I427", "I428", "I429", "I43", "I50", "P290"),
  arrhythmia = c("I441", "I442", "I443", "I456", "I459", "I47", "I48", "I49", 
                 "R000", "R001", "R008", "T821", "Z450", "Z950"),
  valve = c("A520", "I05", "I06", "I07", "I08", "I091", "I098", "I34", "I35", "I36", "I37", "I38", "I39",
            "Q230", "Q231", "Q232", "Q233", "Z952", "Z953", "Z954"),
  pulm.circ = c("I26", "I27", "I280", "I288", "I289"),
  pvd = c("I70", "I71", "I731", "I738", "I739", "I771", "I790", "I792", "K551", "K558", "K559", "Z958", "Z959"),
  htn = c("I10"),
  htncx = c("I11", "I12", "I13", "I15"),
  paralysis = c("G041", "G114", "G801", "G802", "G81", "G82", "G830", "G831", "G832", "G833", "G834", "G839"),
  neuro.other = c("G10", "G11", "G12", "G13", "G20", "G21", "G22", "G254", "G255", "G312", "G318", "G319",
                  "G32", "G35", "G36", "G37", "G40", "G41", "G931", "G934", "R470", "R56"),
  chronic.pulm = c("I278", "I279", "J40", "J41", "J42", "J43", "J44", "J45", "J46", "J47", "J60", "J61", "J62",
                   "J63", "J64", "J65", "J66", "J67", "J684", "J701", "J703"),
  dm.uncomp = c("E100", "E101", "E109", "E110", "E111", "E119", "E120", "E121", "E129", "E130", "E131", "E139",
                "E140", "E141", "E149"),
  dm.comp = c("E102", "E103", "E104", "E105", "E106", "E107", "E108",
              "E112", "E113", "E114", "E115", "E116", "E117", "E118",
              "E122", "E123", "E124", "E125", "E126", "E127", "E128",
              "E132", "E133", "E134", "E135", "E136", "E137", "E138",
              "E142", "E143", "E144", "E145", "E146", "E147", "E148"),
  hypothyroid = c("E00", "E01", "E02", "E03", "E890"),
  renal = c("I120", "I131", "N18", "N19", "N250", "Z490", "Z491", "Z492", "Z940", "Z992"),
  liver = c("B18", "I85", "I864", "I982", "K70", "K711", "K713", "K714", "K715", "K717", "K72", "K73", "K74",
            "K760", "K762", "K763", "K764", "K765", "K766", "K767", "K768", "K769", "Z944"),
  pud = c("K257", "K259", "K267", "K269", "K277", "K279", "K287", "K289"),
  hiv = c("B20", "B21", "B22", "B24"),
  lymphoma = c("C81", "C82", "C83", "C84", "C85", "C88", "C96", "C900", "C902"),
  mets = c("C77", "C78", "C79", "C80"),
  solid.tumor = c("C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12",
                  "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25",
                  "C26", "C30", "C31", "C32", "C33", "C34", "C37", "C38", "C39", "C40", "C41", "C43", "C45",
                  "C46", "C47", "C48", "C49", "C50", "C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58",
                  "C60", "C61", "C62", "C63", "C64", "C65", "C66", "C67", "C68", "C69", "C70", "C71", "C72",
                  "C73", "C74", "C75", "C76", "C97"),
  rheum = c("L940", "L941", "L943", "M05", "M06", "M08", "M120", "M123", "M30", "M310", "M311", "M312", "M313",
            "M32", "M33", "M34", "M35", "M45", "M461", "M468", "M469"),
  coag = c("D65", "D66", "D67", "D68", "D691", "D693", "D694", "D695", "D696"),
  obesity = c("E66"),
  wt.loss = c("E40", "E41", "E42", "E43", "E44", "E45", "E46", "R634", "R64"),
  lytes = c("E222", "E86", "E87"),
  anemia.loss = c("D500"),
  anemia.def = c("D508", "D509", "D51", "D52", "D53"),
  etoh = c("F10", "E52", "G621", "I426", "K292", "K700", "K703", "K709", "T51", "Z502", "Z714", "Z721"),
  drugs = c("F11", "F12", "F13", "F14", "F15", "F16", "F18", "F19", "Z715", "Z722"),
  psychoses = c("F20", "F22", "F23", "F24", "F25", "F28", "F29", "F302", "F312", "F315"),
  depression = c("F204", "F313", "F314", "F315", "F32", "F33", "F341", "F412", "F432")
)

names(icd9_elix_map) <- names(icd10_elix_map) <-
  c("CHF", "Arrhythmia", "Valvular", "PHTN", "PVD", "HTN", "HTNcx", "Paralysis", 
    "NeuroOther", "Pulmonary", "DM", "DMcx", "Hypothyroid", "Renal", 
    "Liver", "PUD", "HIV", "Lymphoma", "Mets", "Tumor", "Rheumatic", 
    "Coagulopathy", "Obesity", "WeightLoss", "FluidsLytes", "BloodLoss", 
    "Anemia", "Alcohol", "Drugs", "Psychoses", "Depression")


icd9_generate_map_elix <- function(save_pkg_data = TRUE) {
  icd9_map_elix <- list(
    chf = c(
      "398.91", "402.11", "402.91", "404.11", "404.13", "404.91",
      "404.93", "428.0" %i9da% "428.9"
    ),
    arrhythmia = c(
      "426.1", "426.11", "426.13", "426.2" %i9da% "426.53",
      "426.6" %i9da% "426.89", "427.0", "427.2", "427.31",
      "427.60", "427.9", "785", "V45.0", "V53.3"
    ),
    valve = c(
      "93.20" %i9da% "93.24", "394.0" %i9da% "397.1",
      "424.0" %i9da% "424.91", "746.3" %i9da% "746.6",
      "V42.2", "V43.3"
    ),
    pulm.circ = c("416.0" %i9da% "416.9", " 417.9"),
    pvd = c(
      "440.0" %i9da% "440.9", "441.2", "441.4", "441.7", "441.9",
      "443.1" %i9da% "443.9", "447.1", "557.1", "557.9", "V43.4"
    ),
    htn = c("401.1", "401.9"),
    htncx = c(
      "402.10", "402.90", "404.10", "404.90", "405.11", "405.19",
      "405.91", "405.99"
    ),
    paralysis = c("342.0" %i9da% "342.12", "342.9" %i9da% "344.9"),
    neuro.other = c(
      "331.9", "332.0", "333.4", "333.5", "334.0" %i9da% "335.9",
      "340", "341.1" %i9da% "341.9", "345.00" %i9da% "345.11",
      "345.40" %i9da% "345.51", "345.80" %i9da% "345.91", "348.1",
      "348.3", "780.3", "784.3"
    ),
    chronic.pulm = c(
      "490" %i9da% "492.8", "493.00" %i9da% "493.91", "494",
      "495.0" %i9da% "505", "506.4"
    ),
    dm.uncomp = c("250.00" %i9da% "250.33"),
    dm.comp = c("250.40" %i9da% "250.73", "250.90" %i9da% "250.93"),
    hypothyroid = c("243" %i9da% "244.2", "244.8", "244.9"),
    renal = c(
      "403.11", "403.91", "404.12", "404.92", "585", "586", "V42.0",
      "V45.1", "V56.0", "V56.8"
    ),
    liver = c(
      "70.32", "70.33", "70.54", "456.0", "456.1", "456.20", "456.21",
      "571.0", "571.2", "571.3", "571.40" %i9da% "571.49", "571.5",
      "571.6", "571.8", "571.9", "572.3", "572.8", "V42.7"
    ),
    pud = c(
      "531.70", "531.90", "532.70", "532.90", "533.70", "533.90",
      "534.70", "534.90", "V12.71"
    ),
    hiv = c("42" %i9da% "44.9"),
    lymphoma = c(
      "200.00" %i9da% "202.38", "202.50" %i9da% "203.01",
      "203.8" %i9da% "203.81", "238.6", "273.3", "V10.71", "V10.72",
      "V10.79"
    ),
    mets = c("196.0" %i9da% "199.1"),
    solid.tumor = c(
      "140.0" %i9da% "172.9", "174.0" %i9da% "175.9",
      "179" %i9da% "195.8", "V10.00" %i9da% "V10.9"
    ),
    rheum = c(
      "701.0", "710.0" %i9da% "710.9", "714.0" %i9da% "714.9",
      "720.0" %i9da% "720.9", "725"
    ),
    coag = c("286.0" %i9da% "286.9", "287.1", "287.3" %i9da% "287.5"),
    obesity = c("278.0"),
    wt.loss = c("260" %i9da% "263.9"),
    lytes = c("276.0" %i9da% "276.9"),
    anemia.loss = c("280.0"),
    anemia.def = c("280.1" %i9da% "281.9", "285.9"),
    etoh = c(
      "291.1", "291.2", "291.5", "291.8", "291.9",
      "303.90" %i9da% "303.93", "305.00" %i9da% "305.03", "V11.3"
    ),
    drugs = c(
      "292.0", "292.82" %i9da% "292.89", "292.9",
      "304.00" %i9da% "304.93", "305.20" %i9da% "305.93"
    ),
    psychoses = c("295.00" %i9da% "298.9", "299.10" %i9da% "299.11"),
    depression = c("300.4", "301.12", "309.0", "309.1", "311")
  )
  icd9_map_elix <- lapply(icd9_map_elix, decimal_to_short.icd9)
  icd9_map_elix <- lapply(
    icd9_map_elix,
    children.icd9,
    short_code = TRUE, defined = FALSE
  )
  names(icd9_map_elix) <- icd::names_elix_htn_abbrev
  icd9_map_elix <- as.comorbidity_map(icd9_map_elix)
  if (save_pkg_data) {
    .save_in_data_dir(icd9_map_elix)
  }
  invisible(icd9_map_elix)
}

# clone the repo (or download zip and unzip)
# git clone https://github.com/jackwasey/icd.git

# then in R:
icd_r_files <- list.files("~/ICU/icd/R", pattern = "\\.R$", full.names = TRUE)
for (f in icd_r_files) source(f, local = TRUE)

