#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_isserlis(IntegerVector ind, NumericVector mu, NumericMatrix Sigma0) {
  int n = ind.size();
  int i, j, k, l;
  
  if (n == 1) {
    i = ind(0);
    return mu(i);
  }
  
  if (n == 2) {
    i = ind(0);
    j = ind(1);
    return mu(i) * mu(j) + Sigma0(i, j);
  }
  
  if (n == 3) {
    i = ind(0);
    j = ind(1);
    k = ind(2);
    return mu(i)*mu(j)*mu(k) +
      mu(i)*Sigma0(j, k) +
      mu(j)*Sigma0(i, k) +
      mu(k)*Sigma0(i, j);
  }
  
  if (n == 4) {
    i = ind(0);
    j = ind(1);
    k = ind(2);
    l = ind(3);
    
    double m4 = mu(i)*mu(j)*mu(k)*mu(l);
    
    double m2s =
      mu(i)*mu(j)*Sigma0(k, l) +
      mu(i)*mu(k)*Sigma0(j, l) +
      mu(i)*mu(l)*Sigma0(j, k) +
      mu(j)*mu(k)*Sigma0(i, l) +
      mu(j)*mu(l)*Sigma0(i, k) +
      mu(k)*mu(l)*Sigma0(i, j);
    
    double covcov =
      Sigma0(i, j)*Sigma0(k, l) +
      Sigma0(i, k)*Sigma0(j, l) +
      Sigma0(i, l)*Sigma0(j, k);
    
    return m4 + m2s + covcov;
  }
  
  // if ind has length > 4 or is empty
  return 1;
}

// [[Rcpp::export]]
double cpp_eta(int i, int j, int k, int l, NumericVector mu, NumericMatrix Sigma0, int p) {
  
  IntegerVector ind;
  double fact = 1;
  
  if (i >= p) ind.push_back(i-p);
  if (j >= p) ind.push_back(j-p);
  if (k >= p) ind.push_back(k-p);
  if (l >= p) ind.push_back(l-p);
  
  if (i >= p && j >= p) fact *= -0.5;
  if (k >= p && l >= p) fact *= -0.5;
  
  if (i >= p && j >= p && i != j) fact *= 2;
  if (k >= p && l >= p && k != l) fact *= 2;

  return fact * cpp_isserlis(ind, mu, Sigma0);
}

// [[Rcpp::export]]
double cpp_lambda(int i, int j, int k, int l, int dim, int p) {
  
  NumericVector ind;
  double ret = 1;
  if (i != -1 && i < p) ret *= ((1<<i & dim) > 0);
  if (j != -1 && j < p) ret *= ((1<<j & dim) > 0);
  if (k != -1 && k < p) ret *= ((1<<k & dim) > 0);
  if (l != -1 && l < p) ret *= ((1<<l & dim) > 0);
  
  if (i != -1 && j != -1 && i < p && j < p && i != j) ret *= 2;
  if (k != -1 && l != -1 && k < p && k < p && k != l) ret *= 2;
  
  return ret;
}

// [[Rcpp::export]]
NumericVector cpp_logc_ij_cts(double Amut, NumericVector Ez, IntegerVector idx, 
                              NumericMatrix Sigma0, NumericMatrix Sig0Lam, 
                              int p, int d) {
  
  NumericVector grad(idx.length());
  NumericVector mu(d);
  double eij;
  int i, j;
  
  for (int x = 0; x < idx.length(); ++x) {
    
    i = (idx(x) / (p+d));
    j = idx(x) % (p+d);
    
    eij = 0;
    
    for (int dim = 0; dim < Ez.length(); ++dim) {
      
      // compute mu
      for (int iz = 0; iz < d; ++iz) {
        mu(iz) = 0;
        for (int ii = 0; ii < p; ++ii) {
          
          mu(iz) += Sig0Lam(iz, ii) * ((1<<ii & dim) > 0);
        }
      }
      
      eij += 1 / Amut * cpp_lambda(i, j, -1, -1, dim, p) * 
        cpp_eta(i, j, -1, -1, mu, Sigma0, p)  * Ez(dim);
    }
    
    grad(x) = eij; 
  }
  
  return grad;
}

// [[Rcpp::export]]
NumericMatrix cpp_hess_cts(double Amut, NumericVector Ez, IntegerVector idx, 
                           NumericMatrix Sigma0, NumericMatrix Sig0Lam, 
                           int p, int d) {
  
  NumericMatrix hess(idx.length(), idx.length());
  NumericVector mu(d);
  double eij, ekl, eijkl;
  int i, j, k, l;
  
  for (int x = 0; x < idx.length(); ++x) {
    
    for (int y = 0; y < idx.length(); ++y) {
      
      i = idx(x) / (p+d);
      j = idx(x) % (p+d);
      k = idx(y) / (p+d);
      l = idx(y) % (p+d);
      
      // Rcpp::Rcout << " x = " << x << " y = " << y << " i = " << i << " j = " << j << " k = " << k << " l = " << l << "\n";
      
      eij = 0, ekl = 0, eijkl = 0;
      
      for (int dim = 0; dim < Ez.length(); ++dim) {
        
        // compute mu
        for (int ai = 0; ai < d; ++ai) {
          mu(ai) = 0;
          for (int bi = 0; bi < p; ++bi) {
            
            mu(ai) += Sig0Lam(ai, bi) * ((1<<bi & dim) > 0);
          }
        }
        
        eij += 1 / Amut * cpp_lambda(i, j, -1, -1, dim, p) *
          cpp_eta(i, j, -1, -1, mu, Sigma0, p) * Ez(dim);
        ekl += 1 / Amut * cpp_lambda(k, l, -1, -1, dim, p) *
          cpp_eta(k, l, -1, -1, mu, Sigma0, p) * Ez(dim);
        eijkl += 1 / (Amut) * cpp_lambda(i, j, k, l, dim, p) *
          cpp_eta(i, j, k, l, mu, Sigma0, p) * Ez(dim);
      }
      
      hess(x, y) = eijkl - eij * ekl;
    }
  }
  
  // returns the Hessian of cumulant H(A) = -H(L), negative Hessian of log-like
  return hess;
}

// [[Rcpp::export]]
NumericVector cpp_pzw_prop(NumericMatrix Sigma, NumericVector LamTw) {
  
  int k = Sigma.ncol(); int pz_dim = 1;
  for (int i = 0; i < k; ++i) pz_dim *= 2;
  NumericVector pz(pz_dim);
  
  for (int i = 0; i < pz_dim; ++i) {
    
    for (int ki = 0; ki < k; ++ki) {
      
      for (int kj = 0; kj < k; ++kj) {
        
        pz(i) += Sigma(ki, kj) * ((i & (1<<ki)) > 0) * ((i & (1<<kj)) > 0) +
          LamTw(ki) * ((i & (1<<kj)) > 0);
      }
    }
    
    pz(i) = exp(pz(i));
  }
  
  return pz;
}

/* ---------- helpers ---------- */
inline double eta_fast(int i,int j,int k,int l,
                       const NumericVector& mu,
                       const NumericMatrix& S,
                       int p){
  /* collect latent indices on the stack */
  int idx[4], n = 0;
  if(i>=p) idx[n++] = i-p;
  if(j>=p) idx[n++] = j-p;
  if(k>=p) idx[n++] = k-p;
  if(l>=p) idx[n++] = l-p;
  
  double f = 1.0;
  if(i>=p && j>=p) f *= -0.5;
  if(k>=p && l>=p) f *= -0.5;
  if(i>=p && j>=p && i!=j) f *= 2.0;
  if(k>=p && l>=p && k!=l) f *= 2.0;
  
  /* Isserlis moments without heap allocations */
  double e = 1.0;
  switch(n){
  case 0: break;
  case 1: e = mu[idx[0]]; break;
  case 2:
    e = mu[idx[0]]*mu[idx[1]] + S(idx[0],idx[1]); break;
  case 3:{
      int a=idx[0],b=idx[1],c=idx[2];
      e = mu[a]*mu[b]*mu[c] +
        mu[a]*S(b,c) + mu[b]*S(a,c) + mu[c]*S(a,b); }
    break;
  case 4:{
    int a=idx[0],b=idx[1],c=idx[2],d=idx[3];
    double m4 = mu[a]*mu[b]*mu[c]*mu[d];
    double m2s = mu[a]*mu[b]*S(c,d) + mu[a]*mu[c]*S(b,d) +
      mu[a]*mu[d]*S(b,c) + mu[b]*mu[c]*S(a,d) +
      mu[b]*mu[d]*S(a,c) + mu[c]*mu[d]*S(a,b);
    double cc  = S(a,b)*S(c,d) + S(a,c)*S(b,d) + S(a,d)*S(b,c);
    e = m4 + m2s + cc; }
    break;
  }
  return f*e;
}

/* ---------- helpers ---------- */
inline double lambda2_fast(int a,int b,int dim,int p,const double* bits){
  double v = 1.0;
  if (a != -1 && a < p){ v *= bits[dim*p + a]; if(v==0.0) return 0.0; }
  if (b != -1 && b < p){ v *= bits[dim*p + b]; if(v==0.0) return 0.0; }
  if (a != -1 && b != -1 && a < p && b < p && a != b) v *= 2.0;
  return v;
}

/* four‑index λ (identical to original logic, incl. the k–l typo) */
inline double lambda4_fast(int i,int j,int k,int l,int dim,int p,const double* bits){
  double v = 1.0;
  if (i != -1 && i < p){ v *= bits[dim*p + i]; if(v==0.0) return 0.0; }
  if (j != -1 && j < p){ v *= bits[dim*p + j]; if(v==0.0) return 0.0; }
  if (k != -1 && k < p){ v *= bits[dim*p + k]; if(v==0.0) return 0.0; }
  if (l != -1 && l < p){ v *= bits[dim*p + l]; if(v==0.0) return 0.0; }
  
  if (i != -1 && j != -1 && i < p && j < p && i != j) v *= 2.0;
  if (k != -1 && l != -1 && k < p && k != l)          v *= 2.0;  // matches old bug
  
  return v;
}


/* ---------- main ---------- */
// [[Rcpp::export]]
NumericMatrix cpp_hess_cts_fast(double Amut, NumericVector Ez, IntegerVector idx,
                                NumericMatrix Sigma0, NumericMatrix Sig0Lam,
                                int p, int d) {
  
  const int Ez_len  = Ez.size();
  const int idx_len = idx.size();
  
  /* μ(a, dim) pre‑compute */
  NumericMatrix mu_mat(d, Ez_len);
  for (int dim = 0; dim < Ez_len; ++dim){
    for (int a = 0; a < d; ++a){
      double s = 0.0;
      for (int b = 0; b < p; ++b)
        if ((dim >> b) & 1) s += Sig0Lam(a, b);
        mu_mat(a, dim) = s;
    }
  }
  
  /* bit‑indicator buffer */
  std::vector<double> bits(Ez_len * p);
  for (int dim = 0; dim < Ez_len; ++dim)
    for (int i = 0; i < p; ++i)
      bits[dim*p + i] = (double)((dim >> i) & 1);
  
  /* Hessian */
  NumericMatrix hess(idx_len, idx_len);
  
  for (int x = 0; x < idx_len; ++x){
    
    int i = idx[x] / (p + d);
    int j = idx[x] % (p + d);
    
    for (int y = 0; y < idx_len; ++y){
      int k = idx[y] / (p + d);
      int l = idx[y] % (p + d);
      
      double eij = 0.0, ekl = 0.0, eijkl = 0.0;
      
      for (int dim = 0; dim < Ez_len; ++dim){
        const NumericVector mu = mu_mat(_, dim);        // view, no copy
        
        double lam_ij   = lambda2_fast(i,j,dim,p,bits.data());
        double lam_kl   = lambda2_fast(k,l,dim,p,bits.data());
        double lam_ijkl = lambda4_fast(i,j,k,l,dim,p,bits.data());
        
        if(lam_ij)   eij   += lam_ij   * eta_fast(i,j,-1,-1,mu,Sigma0,p) * Ez[dim] / Amut;
        if(lam_kl)   ekl   += lam_kl   * eta_fast(k,l,-1,-1,mu,Sigma0,p) * Ez[dim] / Amut;
        if(lam_ijkl) eijkl += lam_ijkl * eta_fast(i,j,k,l,   mu,Sigma0,p) * Ez[dim] / Amut;
      }
      
      hess(x, y) = eijkl - eij * ekl;
    }
  }
  
  return hess;
}
