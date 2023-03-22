
for (x in c(0, 1)) {
  
  for (y in c(0, 1)) {
    
    print(
      colMeans(aux$Z[X == x & Y == y, ]) - 
        colMeans(z_1d(k = 2, 10^4, lambda = lambda[[1+x]][[1+y]]))
    )
  }
}

for (x in c(0, 1)) {
  
  for (y in c(0, 1)) {
    
    print(
      colMeans(aux$Z[X == x & Y == y, ]) - expit(lambda[[1+x]][[1+y]])
    )
  }
}


par(mfrow = c(2, 2))
for (x in c(0, 1)) {
  
  for (y in c(0, 1)) {
    
    dist <- t(do.call(rbind, em[[1+x]][[1+y]]$lambda)) - aux$lxy[[1+x]][[1+y]]
    plot(sqrt(colSums(dist^2)), pch = 19)
  }
}

