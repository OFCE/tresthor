#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix computeJacobian(Rcpp::Function f, NumericVector x, double h ) ;

// Newton-Raphson method for a system of non-linear equations
// [[Rcpp::export]]
NumericVector newtonRaphsonSystem( Rcpp::Function f, NumericVector initialGuess, double tolerance, int maxIter) {
  NumericVector x = clone(initialGuess);
  NumericVector x_new(2);
  NumericVector f_val;
  NumericMatrix jac;
  int iter = 0;

  while (iter < maxIter) {
    f_val = Rcpp::as<NumericVector>(f(x));
    jac = computeJacobian(x);

    // Solve for delta using LU decomposition
    NumericVector delta = solve(jac, f_val);

    // Update the guess
    x_new = x - delta;

    // Check for convergence
    if (sqrt(sum(pow(x_new - x, 2))) < tolerance) {
      return x_new;
    }

    x = x_new;
    iter++;
  }

  Rcpp::warning("Maximum iterations reached without convergence.");
  return x_new;
}
