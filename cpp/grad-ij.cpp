#include <Rcpp.h>
// #include <bits/stdc++.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector cpp_idx(int x, int a, int b, int dim) {
  
  LogicalVector res(1<<dim);
  bool skip;
  x -= 1;
  
  for (int i = 0; i < 1<<dim; ++i) {
    
    skip = false;
    for(int j = 0; j < dim; ++j) {
      
      if ( !(i & (1<<j)) && ((x & (1<<j))))  skip=true;
      // if (i == 1 && skip) {
      //   Rcout<< (1<<j)<<std::endl;
      //   Rcout<< (i & 1<<j)<<std::endl;
      //   Rcout<< (x & 1<<j)<<std::endl;
      // }
    }
    
    // Rcout<< skip<<std::endl;
    if (!skip) {
      res(i) = (i & 1<<(a-1)) && (i & 1<<(b-1));
    }
    
  }
  
  return res;
}

// [[Rcpp::export]]
NumericVector cpp_fi_hash(int x, double fi, int dim) {
  
  NumericVector res(1<<dim);
  bool skip;
  x -= 1;
  
  for (int i = 0; i < 1<<dim; ++i) {
    
    skip = false;
    for(int j = 0; j < dim; ++j) {
      
      if ( !(i & (1<<j)) && ((x & (1<<j)))) {
        res(i) = 0;
        skip = true;
        break;
      }
    }
    
    // Rcout<< skip<<std::endl;
    if (!skip) {
      res(i) = 1;
      for(int j = 0; j < dim; ++j) {
        if ( i & (1<<j) && !(x & (1<<j))) res(i) *= fi;
        if ( x & (1<<j)) res(i) *= 1-fi;
      }
    }
    
  }
  
  return res;
}

// [[Rcpp::export]]
NumericVector cpp_bit_fi_hash(IntegerVector x, double fi) {
  
  int dim = x.length();
  NumericVector res(1<<dim);
  bool skip;
  
  for (int i = 0; i < 1<<dim; ++i) {
    
    skip = false;
    for(int j = 0; j < dim; ++j) {
      
      if ( !(i & (1<<j)) && x(j) ) {
        res(i) = 0;
        skip = true;
        break;
      }
    }
    
    // Rcout<< skip<<std::endl;
    if (!skip) {
      res(i) = 1;
      for(int j = 0; j < dim; ++j) {
        if ( i & (1<<j) && !x(j) ) res(i) *= fi;
        if ( x(j) ) res(i) *= 1-fi;
      }
    }
    
  }

  return res;
}

// [[Rcpp::export]]
NumericVector cpp_grad_ij(NumericVector rhash, int i, int j, int k, double fi,
                          NumericVector pz) {
  
  LogicalVector idx(pz.length());
  NumericVector fi_wgh_hash(pz.length());
  NumericVector grad_ij(rhash.length());

  double idx_sum; double tot_sum;

  for (int x = 0; x < rhash.length(); ++x) {

    fi_wgh_hash = cpp_fi_hash(rhash(x), fi, k);
    idx = cpp_idx(rhash(x), i, j, k);

    idx_sum = 0; tot_sum = 0;
    for (int y = 0; y < idx.length(); ++y) {

      if (idx(y)) idx_sum += fi_wgh_hash(y) * pz(y);
      tot_sum += fi_wgh_hash(y) * pz(y);
    }

    if (tot_sum > 0) grad_ij(x) = idx_sum / tot_sum;
  }

  return grad_ij;
}

// [[Rcpp::export]]
NumericVector cpp_scaling_2d(NumericMatrix Sigma) {
  
  int k = Sigma.ncol(); int pz_dim = 1;
  for (int i = 0; i < k; ++i) pz_dim *= 2;
  NumericVector pz(pz_dim);
  
  for (int i = 0; i < pz_dim; ++i) {
    
    for (int ki = 0; ki < k; ++ki) {
      
      for (int kj = 0; kj < k; ++kj) {
        
        pz(i) += Sigma(ki, kj) * ((i & (1<<ki)) > 0) * ((i & (1<<kj)) > 0);
      }
    }
    
    pz(i) = exp(pz(i));
  }
  
  return pz;
}

// [[Rcpp::export]]
NumericMatrix cpp_hessian_2d_degenerate(NumericVector pz, int p) {
  
  NumericMatrix hess(p * p, p * p);
  double eij, ekl, eijkl;
  int i, j, k, l;
  
  for (int x = 0; x < p*p; ++x) {

    for (int y = 0; y < p*p; ++y) {

      i = x / p;
      j = x % p;
      k = y / p;
      l = y % p;
      // std::cout<<i<<" "<<j<<" "<<k<<" "<<l<<std::endl;
      eij = 0, ekl = 0, eijkl = 0;
      
      for (int dim = 0; dim < pz.length(); ++dim) {
        
        // std::cout<<(1<<i & dim) * ((1<<j & dim) > 0) * pz(dim)<<std::endl;
        eij += ((1<<i & dim) > 0) * ((1<<j & dim) > 0) * pz(dim);
        ekl += ((1<<k & dim) > 0) * ((1<<l & dim) > 0) * pz(dim);
        eijkl += ((1<<i & dim) > 0) * ((1<<j & dim) > 0) * ((1<<k & dim) > 0) * ((1<<l & dim) > 0) * pz(dim);
      }

      hess(x, y) = eijkl - eij * ekl;
      // std::cout<<x<<" "<<y<<std::endl;
    }
  }
  
  return hess;
}

// [[Rcpp::export]]
NumericMatrix cpp_hessian_2d(NumericVector pz, IntegerVector idx, int p) {
  
  NumericMatrix hess(idx.length(), idx.length());
  double eij, ekl, eijkl;
  int i, j, k, l;
  
  for (int x = 0; x < idx.length(); ++x) {
    
    for (int y = 0; y < idx.length(); ++y) {
      
      i = (idx(x) / p);
      j = idx(x) % p;
      k = idx(y) / p;
      l = idx(y) % p;
      // std::cout<<i<<" "<<j<<" "<<k<<" "<<l<<std::endl;
      eij = 0, ekl = 0, eijkl = 0;
      
      for (int dim = 0; dim < pz.length(); ++dim) {
        
        // std::cout<<(1<<i & dim) * ((1<<j & dim) > 0) * pz(dim)<<std::endl;
        eij += ((1<<i & dim) > 0) * ((1<<j & dim) > 0) * pz(dim);
        ekl += ((1<<k & dim) > 0) * ((1<<l & dim) > 0) * pz(dim);
        eijkl += ((1<<i & dim) > 0) * ((1<<j & dim) > 0) * ((1<<k & dim) > 0) * ((1<<l & dim) > 0) * pz(dim);
      }
      
      // if (x == 0 & y == 0) {
      //   std::cout<<i<<" "<<j<<std::endl;
      //   // std::cout<<(1 + i!=j)<<" "<<(1 + k!=l)<<std::endl;
      //   std::cout<<eij<<" "<<ekl<<" "<<eijkl<<std::endl;
      // }
      hess(x, y) = (eijkl - eij * ekl) * (1 + (i!=j)) * (1 + (k!=l));
      // std::cout<<x<<" "<<y<<std::endl;
    }
  }
  
  return hess;
}

// [[Rcpp::export]]
NumericMatrix cpp_mu(NumericVector pz, int p) {
  
  NumericMatrix mu(p, p);
  
  for (int i = 0; i < p; ++i) {

    for (int j = 0; j < p; ++j) {
      
      for (int dim = 0; dim < pz.length(); ++dim) {
        
        // std::cout<<(1<<i & dim) * ((1<<j & dim) > 0) * pz(dim)<<std::endl;
        mu(i, j) += ((1<<i & dim) > 0) * ((1<<j & dim) > 0) * pz(dim);
      }
    }
  }
  
  return mu;
}