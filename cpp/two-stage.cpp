#include <Rcpp.h>
using namespace Rcpp;

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