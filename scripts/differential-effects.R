
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))


src <- c("miiv", "aumc", "hirid", "anzics")
plt <- list()
for (dsrc in src) {
  
  dat <- load_data(dsrc)
  mod_uni <- fit_log(dat, "univar")
  mod_mul <- fit_log(dat, "multivar")
  
  plt[[dsrc]] <- cowplot::plot_grid(
    plot_or(mod_uni) + ggtitle(paste(srcwrap(dsrc), "univariate")), 
    plot_or(mod_mul) + ggtitle(paste(srcwrap(dsrc), "multivariate"))
  )
}

cowplot::plot_grid(plotlist = plt, nrow = length(plt))

ggsave(file.path(root, "results", "sex-impact-with-age.png"), width = 10, height = 14)

#' * differential impact of sex *
#
# tags <- c("<18", "<25", "<30", "<35", ">35")
# ggplot(dat[, list(mean(death), .N), by = c("sex", "bmi_bin")],
#        aes(x = bmi_bin, y = V1, color = sex, fill = sex)) +
#   geom_point() + geom_line() +
#   geom_ribbon(aes(ymin = V1 - 1.96 * V1 * (1-V1) / sqrt(N),
#                   ymax = V1 + 1.96 * V1 * (1-V1) / sqrt(N)), alpha = 0.4) +
#   scale_fill_discrete(guide = "none") +
#   scale_x_continuous(labels = tags, breaks = seq_along(tags)) +
#   theme_bw() + xlab("BMI group") + ylab("Mortality rate") +
#   scale_y_continuous(labels = scales::percent) +
#   theme(legend.position = c(0.8, 0.8),
#         legend.box.background = element_rect())
# 
# ggsave(file.path(root, "results", "sex.png"), width = 5, height = 5)

#' * differential impact of height *
# 
# med_height <- median(dat$height)
# med_height_m <- median(dat[sex == "Male"]$height)
# dat[, height_bin := height > med_height]
# dat[, height_bin_m := height > med_height_m]
# 
# height_all <- ggplot(dat[, list(mean(death), .N), by = c("height_bin", "bmi_bin")],
#                      aes(x = bmi_bin, y = V1, color = height_bin, fill = height_bin)) +
#   geom_point() + geom_line() +
#   geom_ribbon(aes(ymin = V1 - 1.96 * V1 * (1-V1) / sqrt(N),
#                     ymax = V1 + 1.96 * V1 * (1-V1) / sqrt(N)), alpha = 0.4) +
#   scale_fill_discrete(guide = "none") +
#   scale_x_continuous(labels = tags, breaks = seq_along(tags)) +
#   theme_bw() + xlab("BMI group") + ylab("Mortality rate") +
#   scale_y_continuous(labels = scales::percent) +
#   scale_color_discrete(name = "Height", labels = c("Below Median", "Above Median")) +
#   theme(legend.position = c(0.8, 0.8),
#         legend.box.background = element_rect()) + ggtitle("All patients")
# 
# height_male <- ggplot(dat[sex == "Male", list(mean(death), .N), by = c("height_bin_m", "bmi_bin")],
#                       aes(x = bmi_bin, y = V1, color = height_bin_m, fill = height_bin_m)) +
#   geom_point() + geom_line() +
#   geom_ribbon(aes(ymin = V1 - 1.96 * V1 * (1-V1) / sqrt(N),
#                   ymax = V1 + 1.96 * V1 * (1-V1) / sqrt(N)), alpha = 0.4) +
#   scale_fill_discrete(guide = "none") +
#   scale_x_continuous(labels = tags, breaks = seq_along(tags)) +
#   theme_bw() + xlab("BMI group") + ylab("Mortality rate") +
#   scale_y_continuous(labels = scales::percent) +
#   scale_color_discrete(name = "Height", labels = c("Below Median", "Above Median")) +
#   theme(legend.position = c(0.8, 0.8),
#         legend.box.background = element_rect()) + ggtitle("Male only")
# 
# cowplot::plot_grid(height_all, height_male, ncol = 2L)
# 
# ggsave(file.path(root, "results", "height.png"), width = 10, height = 5)

#' * within BMI comorbidities *

# dat[, with_norm := .bincode(bmi, c(18, 21, 25))]
# 
# chrl <- dat[, rowSums(.SD), .SDcols = config("charlson")[[1]]]
# 
# dat$charlson <- chrl
# 
# boot <- NULL
# for (i in 1:100) {
#   idx <- sample(nrow(dat), replace = TRUE)
#   boot <- rbind(boot, dat[idx, list(chrl = mean(charlson), rep = i), by = c("bmi_bin")])
# }
# 
# charlson_bin <- ggplot(boot[, list(mean(chrl), sd(chrl)), by = c("bmi_bin")],
#        aes(x = bmi_bin, y = V1)) +
#   geom_ribbon(aes(ymin = V1 - 1.96 * V2, ymax = V1 + 1.96 * V2), alpha = 0.4) +
#   geom_point() + geom_line() +
#   scale_x_continuous(labels = tags, breaks = seq_along(tags)) +
#   theme_bw() + xlab("BMI group") + ylab("# comorbidities (Charlson)") #+
  # scale_y_continuous(labels = scales::percent) #+
  # scale_color_discrete(name = "Height", labels = c("Below Median", "Above Median")) +
  # theme(legend.position = c(0.8, 0.8),
  #       legend.box.background = element_rect())


#' * a LOESS smoother *
# charlson_loess <- ggplot(dat[bmi < 40 & bmi> 15], aes(x = bmi, y = charlson)) +
#   geom_smooth() + xlab("BMI value") + ylab("# comorbidities (Charlson)") +
#   theme_bw()
# 
# cowplot::plot_grid(charlson_bin, charlson_loess, ncol = 2L)
# 
# ggsave(file.path(root, "results", "charlson.png"), width = 10, height = 5)
