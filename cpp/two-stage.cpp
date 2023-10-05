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