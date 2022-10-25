
evl_frm <- function(fwd, aux) {
  
  data.frame(adist = unlist(fwd$adverse), 
             bdist = abs(unlist(fwd$beta) - aux$beta),
             zdist = unlist(
               lapply(fwd$pz, function(z) 100 * sum(abs(z - aux$pz)))
             ),
             iter = seq_along(fwd$adverse))
}

vis_evl <- function(evl) {
  
  ggplot(reshape2::melt(evl, id.vars = "iter"), aes(x = iter, y = value)) +
    geom_line() + geom_point() + theme_bw() +
    facet_grid(rows = vars(variable), scales = "free")
}



