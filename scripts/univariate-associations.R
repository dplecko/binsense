library(ricu)
library(icd)
library(ranger)
library(latex2exp)
library(ggplot2)
library(stringr)
library(cowplot)

root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

src <- "mimic"
dat <- load_concepts(c("death", "diag", "elix", "bmi_bin"), src)
dat[, c(index_var(dat)) := NULL]
dat <- dat[!is.na(bmi_bin)]
dat[is.na(death), "death"] <- FALSE
dat[is.na(diag), "diag"] <- "OTH"

acu <- load_concepts("sofa", src, explicit_wins = hours(24L))
acu <- acu[, c(index_var(acu)) := NULL]
acu <- rename_cols(acu, "acu", "sofa")

dat <- merge(dat, acu, all.x = TRUE)
dat <- rename_cols(dat, c("Death", "SOFA", "Elixhauser", "BMI", "Diagnosis"),
                   c("death", "acu", "elix", "bmi_bin", "diag"))

plt <- list()
for (i in 1+seq_along(names(dat)[-1])) {
  for (j in 1+seq_along(names(dat)[-1])) {
    
    if (i == j) p <- ggplot() + theme_void() else {
      
      
      if (length(unique(dat[[i]])) == 2) {
        fill <- names(dat)[i]
        x <- names(dat)[j]
      } else if (length(unique(dat[[j]])) == 2) {
        fill <- names(dat)[j]
        x <- names(dat)[i]
      }
      
      if (is.numeric(dat[[x]])) {
        p <- ggplot(dat, aes_string(x=x, fill=fill)) +
          geom_histogram(aes(y = stat(density*width)), position='dodge') +
          scale_y_continuous(labels = scales::percent) +
          ylab("Proportion")  +
          theme(
            legend.position = c(0.8, 0.75),
            legend.box.background = element_rect("black")
          )
      } else {
        p <- ggplot(dat, aes_string(x=x, fill=fill)) +
          geom_histogram(aes(y = stat(count)), position='dodge', 
                         stat="count") + ylab("Total number") +
          theme(
            legend.position = c(0.8, 0.75),
            legend.box.background = element_rect("black")
          )
      }

    }
    
    plt <- c(plt, list(p))

  }
}

plot_grid(plotlist = plt, ncol = length(names(dat)[-1]))
ggsave("associations.png", width = 20, height = 15)
