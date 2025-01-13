#include <Rcpp.h>
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
    }
    
    if (!skip) {
      res(i) = (i & 1<<(a-1)) && (i & 1<<(b-1));
    }
    
  }
  
  return res;
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
      
      eij = 0, ekl = 0, eijkl = 0;
      
      for (int dim = 0; dim < pz.length(); ++dim) {
        
        eij += ((1<<i & dim) > 0) * ((1<<j & dim) > 0) * pz(dim);
        ekl += ((1<<k & dim) > 0) * ((1<<l & dim) > 0) * pz(dim);
        eijkl += ((1<<i & dim) > 0) * ((1<<j & dim) > 0) * ((1<<k & dim) > 0) * ((1<<l & dim) > 0) * pz(dim);
      }
      
      hess(x, y) = (eijkl - eij * ekl) * (1 + (i!=j)) * (1 + (k!=l));
    }
  }
  
  return hess;
}

// [[Rcpp::export]]
NumericVector cpp_p_z(NumericVector lambda) {
  
  int k = lambda.length(); int pz_dim = 1;
  for (int i = 0; i < k; ++i) pz_dim *= 2;
  NumericVector px_z(pz_dim);
  
  for (int i = 0; i < pz_dim; ++i) {
    
    for (int ki = 0; ki < k; ++ki) {
      
      px_z(i) += lambda(ki) * ((i & (1<<ki)) > 0);
    }
  }
  
  return px_z;
}

// [[Rcpp::export]]
NumericVector cpp_idx_to_bit(IntegerVector idx, int dim) {
  
  NumericMatrix r_bit(idx.length(), dim);
  
  for (int i = 0; i < idx.length(); ++i) {
    
    for (int ki = 0; ki < dim; ++ki) {
      
      r_bit(i, ki) = (idx(i) & (1<<ki)) > 0;
    }
  }
  
  return r_bit;
}

// [[Rcpp::export]]
NumericVector cpp_bit_to_idx_mat(NumericMatrix bit) {
  
  NumericVector idx(bit.nrow());
  int dim = bit.ncol();
  
  for (int i = 0; i < bit.nrow(); ++i) {
    
    for (int ki = 0; ki < dim; ++ki) {
      
      idx(i) += bit(i, ki) * (1<<ki);
    }
  }
  
  return idx;
}

// [[Rcpp::export]]
int cpp_bit_to_idx(NumericVector bit) {
  
  int idx;
  int dim = bit.length();
  idx = 0;
  
  for (int i = 0; i < bit.length(); ++i) {
    
    idx += bit(i) * (1<<i);
  }
  
  return idx;
}
