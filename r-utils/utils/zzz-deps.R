
pkgs <- c("ricu", "readr", "ranger", "ggplot2", "icd", "plotly",
          "latex2exp", "stringr", "assertthat", "data.table", "scales",
          "faircause", "Rcpp", "parallel", "foreach", "doParallel")

Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config"))
Sys.setenv("RICU_SRC_LOAD" = 
             "mimic,miiv,aumc,hirid,eicu,eicu_demo,mimic_demo,anzics")

# if (!all(vapply(pkgs, requireNamespace, logical(1L)))) {
#   stop("Packages {pkgs} are required in order to proceed.")
#   if (!interactive()) q("no", status = 1, runLast = FALSE)
# }

library(readr)
library(ranger)
library(ggplot2)
library(latex2exp)
library(stringr)
library(grf)
library(assertthat)
library(data.table)
library(plotly)
library(ricu)
library(icd)
library(scales)
library(faircause)
library(Rcpp)
library(parallel)
library(binsensate)
library(foreach)
library(doParallel)
library(future)
library(future.apply)

