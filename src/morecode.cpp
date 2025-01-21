#include <cpp11.hpp>
#include <cpp11eigen.hpp>

using namespace Eigen;
using namespace cpp11;

[[cpp11::register]] // allows using the function in R
doubles_matrix<>
solve_mat(doubles_matrix<> x) {
  MatrixXd Y = as_Matrix(x);      // convert from R to C++
  MatrixXd Yinv = Y.inverse();    // Y^(-1)
  return as_doubles_matrix(Yinv); // convert from C++ to R
}
