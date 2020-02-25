#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP quadFormByRow(const Eigen::MappedSparseMatrix< double > S,
                       const Eigen::MappedSparseMatrix< double > X) {
  return wrap(((X * S).cwiseProduct(X)) * Eigen::VectorXd::Ones(X.cols()));
}

