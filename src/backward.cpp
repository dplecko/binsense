#include <Rcpp.h>
using namespace Rcpp;

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

// [[Rcpp::export]]
NumericVector cpp_fi_grad(NumericMatrix Ainv, 
                          NumericVector neg_idx, NumericVector pr, 
                          NumericVector fi,
                          int verbose = 0) {
  
  int k = fi.length();
  NumericVector fi_grad(k);
  int i_bit, j_bit, i_bitb, j_bitb, bit_diff;
  int i;
  // int flag =  0;
  double der_ij, der_ij_new, der_ij_t2;
  
  for (int ki = 0; ki < k; ++ki) {
    
    // comment loop over neg_idx
    for (int i_idx = 0; i_idx < neg_idx.length(); ++i_idx) {
      
      i = neg_idx(i_idx) - 1;
      // need to convert the i...
      
      // loop over 2^k possible r values
      for (int j = 0; j < (1<<k); ++j) {
        
        bit_diff = 0;
        if (Ainv(i, j) == 0) continue;
        
        // is bit in z^(i)?
        i_bit = (i & (1<<ki)) > 0;
        
        // is bit in r^(j) \ z^(i)?
        j_bit = (j & (1<<ki)) > 0;
        
        if (i_bit) {
          
          der_ij = Ainv(i, j) / (1-fi(ki)); // it was a mistake here!!
          
          der_ij_new = 1;
          
          for (int b = 0; b < k; ++b) {
            
            i_bitb = (i & (1<<b)) > 0;
            j_bitb = (j & (1<<b)) > 0;
            bit_diff += j_bitb - i_bitb;
            
            if (j_bitb == 1 & i_bitb == 0) der_ij_new *= fi(b) / (1- fi(b));
            if (j_bitb == 1 & i_bitb == 1) {
              
              der_ij_new *= 1 / (1 - fi(b));
              if (b == ki) {
                
                der_ij_new *= 1 / (1 - fi(b));
              }
            }
          }
          
          if (bit_diff % 2) der_ij_new *= -1;
        } else if (j_bit) {
          
          der_ij = Ainv(i, j) / fi(ki) + Ainv(i, j) / (1 - fi(ki));
          
          der_ij_new = 1;
          
          for (int b = 0; b < k; ++b) {
            
            i_bitb = (i & (1<<b)) > 0;
            j_bitb = (j & (1<<b)) > 0;
            bit_diff += j_bitb - i_bitb;
            
            if (j_bitb == 1 & i_bitb == 0) der_ij_new *= fi(b) / (1 - fi(b));
            if (j_bitb == 1 & i_bitb == 1) {
              
              der_ij_new *= 1 / (1 - fi(b));
              if (b == ki) {
                
                der_ij_new *= 1 / (1 - fi(b));
              }
            }
          }
          
          der_ij_t2 = 1;
          for (int b = 0; b < k; ++b) {
            
            i_bitb = (i & (1<<b)) > 0;
            j_bitb = (j & (1<<b)) > 0;
            
            if (j_bitb == 1 & i_bitb == 0 & b != ki) der_ij_t2 *= fi(b) / (1 - fi(b));
            if (j_bitb == 1 & i_bitb == 0 & b == ki) der_ij_t2 *= 1 / (1 - fi(b));
            if (j_bitb == 1 & i_bitb == 1 & b != ki) der_ij_t2 *= 1 / (1 - fi(b));
          }
          
          der_ij_new += der_ij_t2;
          
          if (bit_diff % 2) der_ij_new *= -1;
        }
        
        fi_grad(ki) += der_ij * pr(j);
      }
    }
  }
  
  return fi_grad;
}
