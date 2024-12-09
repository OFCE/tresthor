// [[Rcpp::depends("RcppArmadillo")]]
// # include <RcppArmadillo.h>
# include <Rcpp.h>

// using namespace arma;
using namespace Rcpp;
// using namespace std;


// [[Rcpp::export]]
NumericMatrix computeJacobian(Rcpp::Function f, NumericVector x, int time ,  double h = 1e-5) {


  NumericVector fx = Rcpp::as<NumericVector>(f( time)); // Evaluate function at x

  int n = x.size(); // Dimension of input vector
  int m = fx.size(); // Dimension of output vector
  NumericMatrix J(m, n); // Initialize Jacobian matrix

  for (int j = 0; j < n; ++j) {
    // Create a vector for modified input x+h
    NumericVector xh = clone(x);
    xh[j] += h; // Increment x[j] by h
    NumericVector fxh = Rcpp::as<NumericVector>(f(t = time));; // Evaluate function at x+h

    // Compute partial derivatives using finite difference
    for (int i = 0; i < m; ++i) {
      J(i, j) = (fxh[i] - fx[i]) / h;
    }
  }

  return J;
}

// // [[Rcpp::depends("RcppArmadillo")]]
// # include <RcppArmadillo.h>
// # include <Rcpp.h>
//
// using namespace arma;
// using namespace Rcpp;
// using namespace std;
//
// NumericMatrix computeJacobian(Rcpp::Function f, NumericVector x, double h);
//
// // [[Rcpp::export]]
// arma::vec newtonRaphsonSystem( Rcpp::Function f, arma::vec initialGuess, double tolerance, int maxIter, int time) {
//   arma::vec x = clone(initialGuess);
//   arma::vec x_new(2);
//   arma::vec f_val;
//   mat jac;
//   int iter = 0;
//
//   while (iter < maxIter) {
//     f_val = Rcpp::as<NumericVector>(f(x));
//     jac = computeJacobian(f, x);
//
//     // Solve for delta using LU decomposition
//     NumericVector delta = Rcpp::solve(jac, f_val);
//
//     // Update the guess
//     x_new = x - delta;
//
//     // Check for convergence
//     if (sqrt(sum(pow(x_new - x, 2))) < tolerance) {
//       return x_new;
//     }
//
//     x = x_new;
//     iter++;
//   }
//
//   Rcpp::warning("Maximum iterations reached without convergence.");
//   return x_new;
// }
