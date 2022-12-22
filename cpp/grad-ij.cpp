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

      // Rcout<<y<<idx(y)<<fi_wgh_hash(y)<<pz(y)<<std::endl;

      if (idx(y)) idx_sum += fi_wgh_hash(y) * pz(y);
      tot_sum += fi_wgh_hash(y) * pz(y);
    }

    if (tot_sum > 0) grad_ij(x) = idx_sum / tot_sum;
  }

  return grad_ij;
}