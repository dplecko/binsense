#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix cpp_A_xy(double fi_xy, int k) {
  
  NumericMatrix Axy(1<<k, 1<<k);
  int i_bit, j_bit;
  
  for (int i = 0; i < (1<<k); ++i) {
    
    for (int j = 0; j < (1<<k); ++j) {
      
      Axy(i, j) = 1;
      
      for (int b = 0; b < k; ++b) {
        
        i_bit = (i & (1<<b)) > 0;
        j_bit = (j & (1<<b)) > 0;
        
        if (i_bit > j_bit) {
          Axy(i, j) = 0;
          continue;
        }
        
        if (i_bit == 1) Axy(i, j) *= (1 - fi_xy);
        if (j_bit - i_bit == 1) Axy(i, j) *= fi_xy;
      }
      
    }
  }
  
  return Axy;
}

// [[Rcpp::export]]
NumericMatrix cpp_invert_A_xy(NumericMatrix A, NumericVector pz, 
                              NumericVector pr) {
  
  NumericMatrix B(A.nrow(), A.ncol());
  
  for (int i = 0; i < B.nrow(); ++i) {
    
    for (int j = 0; j < B.ncol(); ++j) {
      
      B(i, j) = A(j, i) * pz(i) / pr(j);
    }
  }
  
  return B;
}
