
#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/MatrixFunctions>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


RcppExport SEXP mpow( SEXP Q, SEXP t )
{
BEGIN_RCPP
  Eigen::MatrixXd A = Rcpp::as<Eigen::MatrixXd>(Q);
  double m = Rcpp::as<double>(t);
  return Rcpp::wrap(A.pow(m));
END_RCPP
}

