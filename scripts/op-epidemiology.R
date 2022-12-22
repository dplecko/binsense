
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config"))

src <- "miiv"
pids <- id_col(load_concepts(c("bmi"), "miiv")[bmi > 18.5])
cnc <- c("diag", "death", "elix", "DM", "smoke", "age", "bmi_bin",
         "is_hypo", "tw_avg_glu")
dat <- load_concepts(cnc[1L], src, patient_ids = pids)
for (c in cnc[-1L]) {
  dat <- merge(dat, load_concepts(c, src), all.x = TRUE)
}
# dat <- load_concepts(cnc, src)
dat[, c(index_var(dat)) := NULL]
dat <- dat[!is.na(bmi_bin)]
dat[is.na(diag), "diag"] <- "OTH"
dat[is.na(death), "death"] <- FALSE
dat[is.na(DM), "DM"] <- 0
dat[is.na(smoke), "smoke"] <- 0
dat[is.na(elix), "elix"] <- 0
dat[is.na(is_hypo), "is_hypo"] <- FALSE
dat[is.na(tw_avg_glu), "tw_avg_glu"] <- 108


acu <- load_concepts("sofa", src, explicit_wins = hours(24L))
acu <- acu[, c(index_var(acu)) := NULL]
acu <- rename_cols(acu, "acu", "sofa")

dat <- merge(dat, acu, all.x = TRUE)
#dat <- dat[get(id_var(dat)) %in% id_col(load_concepts("DM", "mimic")[DM == 1])]

res <- CausalExplanation_TV(dat, X = "bmi_bin", W = c("diag", "acu", "is_hypo",
                                                      "tw_avg_glu"),
                            Z = c("elix", "age", "DM", "smoke"), Y = "death", 
                            x0 = "<25", x1 = ">25",
                            probability = TRUE)

TV_family(res, "MIMIC-4 (non-thin)")
ggsave("causal_mediation.png", width = 12, height = 7, bg = "white")
