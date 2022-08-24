#define ARMA_DONT_PRINT_ERRORS
#define _USE_MATH_DEFINES

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

/*
 * Code derived from GPvecchia package (Katzfuss et al.)
 */

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

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
arma::mat rev_matrix(arma::mat x) {
  return reverse(x, 1);
}

double sqdist(rowvec l1, rowvec l2) { 
  double ssq = 0.0;
  for(arma::uword k = 0; k < l1.size(); ++k) {
    ssq += (l1[k] - l2[k])*(l1[k] - l2[k]);
  }
  return ssq;
}

arma::mat calc_sqdist(arma::mat x) {
  arma::uword outrows = x.n_rows;
  arma::uword outcols = x.n_rows;
  arma::mat out(outrows, outcols);
  for (arma::uword arow = 0 ; arow < outrows ; arow++) {
    for (arma::uword acol = 0 ; acol < outcols ; acol++) {
      out(arow, acol) = sqdist(x.row(arow), x.row(acol));
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat Exp2Fun_SW(arma::mat distmat, arma::vec covparms, arma::vec g) { 
  // distmat = matrix of SQUARED distances
  // covparms = c(tau2, theta, g, v = NULL)
  int d1 = distmat.n_rows;
  int d2 = distmat.n_cols;
  int j1, j2;
  arma::mat covmat(d1, d2);
  double r;
  for (j1 = 0; j1 < d1; j1++) {
    for (j2 = 0; j2 < d2; j2++) {
      r = distmat(j1, j2)/covparms(1);
      covmat(j1, j2) = covparms(0)*exp(-r);
    }
  }
  if (d1 == d2) {
    for (j1 = 0; j1 < d1; j1++) 
      covmat(j1, j1) += covparms(0) * g(j1);
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat MaternFun_SW(arma::mat distmat, arma::vec covparms, arma::vec g) { 
  // distmat = matrix of SQUARED distances
  // covparms = c(tau2, theta, g, v)
  int d1 = distmat.n_rows;
  int d2 = distmat.n_cols;
  int j1, j2;
  arma::mat covmat(d1, d2);
  double r;
  if (covparms(3) == 0.5) { 
    for (j1 = 0; j1 < d1; j1++) {
      for (j2 = 0; j2 < d2; j2++) {
        r = sqrt(distmat(j1, j2) / covparms(1));
        covmat(j1, j2) = covparms(0) * exp(-r);
      }
    }
  } else if(covparms(3) == 1.5) {
    for (j1 = 0; j1 < d1; j1++) {
      for (j2 = 0; j2 < d2; j2++) {
        r = sqrt(3 * distmat(j1, j2) / covparms(1));
        covmat(j1, j2) = covparms(0) * (1 + r) * exp(-r);
      }
    }
  } else if(covparms(3) == 2.5) {
    for (j1 = 0; j1 < d1; j1++) {
      for (j2 = 0; j2 < d2; j2++) {
        r = sqrt(5 * distmat(j1, j2) / covparms(1));
        covmat(j1, j2) = covparms(0) * (1 + r + pow(r, 2) / 3) * exp(-r);
      }
    }
  } 
  if (d1 == d2) {
    for (j1 = 0; j1 < d1; j1++) 
      covmat(j1, j1) += covparms(0) * g(j1);
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat U_entries_SW(const int Ncores, const arma::uword n, const arma::mat& locs, 
                     const arma::umat& revNNarray, const arma::mat& revCondOnLatent, 
                     const arma::vec covparms, const arma::vec& g){
  const uword m = revNNarray.n_cols - 1;
  const uword Nlocs = locs.n_rows;
  arma::mat Lentries = zeros(Nlocs, m + 1);
  
  #ifdef _OPENMP
    
  #pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
    
    for (uword k = 0; k < Nlocs; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::vec revCon_row = revCondOnLatent.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat dist = calc_sqdist(locs.rows(inds00));
      arma::vec gsub = g.elem(inds00);
      arma::mat covmat = MaternFun_SW(dist, covparms, gsub);
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(k, span(0, n0 - 1)) = M.t();
    }
    
  #else
    
    for (uword k = 0; k < Nlocs; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::vec revCon_row = revCondOnLatent.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat dist = calc_sqdist(locs.rows(inds00));
      arma::vec gsub = g.elem(inds00);
      arma::mat covmat = MaternFun_SW(dist, covparms, gsub);
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(k, span(0, n0 - 1)) = M.t();
    }
    
  #endif
  
  return Lentries;
}
