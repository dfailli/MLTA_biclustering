#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericVector callLogisticRegressionb(arma::mat yy,
                                           arma::mat xxb,
                                           arma::mat offsb,
                                           arma::mat ww22,
                                           arma::vec coeffb) {
  
  int QGMN = xxb.n_rows;
  int P = xxb.n_cols;
  
  int maxIterations = 100;
  double convergenceThreshold = 1e-6;
  
  arma::mat onesMatrix(QGMN, 1, arma::fill::ones);
  
  for (int iteration = 0; iteration < maxIterations; iteration++) {
    arma::mat linp = xxb * coeffb + offsb;
    arma::mat pi = onesMatrix / (onesMatrix + exp(-linp));
    
    arma::vec gradient = xxb.t() * (ww22 % (yy - pi));
    
    arma::mat hessian(P, P, arma::fill::zeros);
    
    for (int i = 0; i < QGMN; i++) {
      arma::rowvec xx_i = xxb.row(i);
      arma::rowvec pdis = pi.row(i);
      arma::rowvec we = ww22.row(i);
      
      arma::mat Fi = xx_i.t() * (we % (arma::diagmat(pdis) - pdis.t() * pdis)) * xx_i;
      hessian += Fi;
    }
    
    arma::vec delta = solve(hessian, gradient, solve_opts::fast);
    
    coeffb += delta;
    
    double deltaNorm = norm(delta, 2);
    
    if (deltaNorm < convergenceThreshold) {
      break;
    }
  }
  
  return wrap(coeffb);
}