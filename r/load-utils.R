
load_difftime.anzics_tbl <- function(x, rows, cols = colnames(x),
                                     id_hint = id_vars(x),
                                     time_vars = ricu::time_vars(x), ...) {
  
  ricu:::warn_dots(...)

  ricu:::load_mihi(x, {{ rows }}, cols, id_hint, time_vars)
}
