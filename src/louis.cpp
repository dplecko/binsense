#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cpp_hess_lambda(NumericVector pz, NumericVector px_z, int p) {
  
  NumericMatrix hess(p+1, p+1);
  int zx, zy;
  
  for (int x = 0; x < p+1; ++x) {
    
    for (int y = 0; y < p+1; ++y) {
      
      for (int i = 0; i < pz.length(); ++i) {
        
        if (x == 0) {
          zx = 1;
        } else {
          zx = (1<<(x-1) & i) > 0;
        }
        
        if (y == 0) {
          zy = 1;
        } else {
          zy = (1<<(y-1) & i) > 0;
        }
        
        hess(x, y) += pz(i) * px_z(i) * (1 - px_z(i)) * zx * zy;  
      }
    }
  }
  
  return hess;
}

// [[Rcpp::export]]
NumericMatrix cpp_hess_mu(NumericVector pz, NumericVector px_z, 
                          NumericVector py_x0z, NumericVector py_x1z, int p) {
  
  NumericMatrix hess(p+2, p+2);
  int zx, zy;
  
  for (int x = 0; x < p+2; ++x) {
    
    for (int y = 0; y < p+2; ++y) {
      
      for (int i = 0; i < pz.length(); ++i) {
        
        for (int xval = 0; xval < 2; ++xval) {
          
          if (x == 0) {
            zx = 1;
          } else if (x == p+1) {
            zx = xval;
          } else {
            zx = (1<<(x-1) & i) > 0;
          }
          
          if (y == 0) {
            zy = 1;
          } else if (y == p+1) {
            zy = xval;
          } else {
            zy = (1<<(y-1) & i) > 0;
          }
          
          if (xval == 0) {
            
            hess(x, y) += pz(i) * (1 - px_z(i)) * py_x0z(i) * (1 - py_x0z(i)) * zx * zy;
          } else {
            
            hess(x, y) += pz(i) * px_z(i) * py_x1z(i) * (1 - py_x1z(i)) * zx * zy;
          }
        }  
      }
    }
  }
  
  return hess;
}
