
# understand the \gamma for GS in vicinity of a random th_s
th_s <- gen_params(k, rand = TRUE)
cov_evals <- evals_Sigma(th_s$Sigma)

# get a function that perturbs the \theta parameter
res <- NULL
for (size in c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1)) {
  
  for (i in 1:100) {
    
    th_t <- perturb_th(th_s, size = size)
    
    gamma_p <- vec_norm(nabla_Q(th_t, th_t, th_s) - nabla_Q(th_t, th_s, th_s)) /
                  vec_norm(em_params(th_t) - em_params(th_s))
    
    th_rand <- gen_params(k, rand = TRUE)
    th_rand$Sigma <- th_s$Sigma
    gamma_p1 <- vec_norm(nabla_Q(th_rand, th_t, th_s) - nabla_Q(th_rand, th_s, th_s)) / 
      vec_norm(em_params(th_t) - em_params(th_s))
    
    gamma_csb <- csb_bound(th_t, th_t, th_s) / 
      vec_norm(em_params(th_t) - em_params(th_s))
    
    res <- rbind(res, c(size, i, gamma_p, gamma_r, gamma_csb))
  }
}

res <- as.data.frame(res)
names(res) <- c("size", "rep", "P", "RAND", "CSB")

res <- reshape2::melt(res, id.vars = c("size", "rep"))

ggplot(res[res$variable %in% c("P", "RAND", "CSB"),], aes(x = value, fill = variable)) +
  geom_density(alpha = 0.4) + theme_bw() +
  facet_wrap(~ size, nrow = length(unique(res$size))) +
  geom_vline(xintercept = min(cov_evals), linetype = "dashed", color = "red")
