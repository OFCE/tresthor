#include <Rcpp.h>

// Define an Rcpp function
// [[Rcpp::export]]
double square(double x) {
  return x * x;
}

#include <Rcpp.h>

// Declare the square function
double square(double x);

// Define another Rcpp function that uses the square function
// [[Rcpp::export]]
double calculateSquareSum(double a, double b) {
  // Call the square function within the Rcpp function
  double squareA = square(a);
  double squareB = square(b);

  // Sum the squared values
  double result = squareA + squareB;

  return result;
}
