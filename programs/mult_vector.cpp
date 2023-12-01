/*
 * Calculates the product of square matrix and vector.
 */
#include <iostream>
#include "linear.h"
#include <omp.h>

int main(void) {
  int dim;
  double start, end;

  std::cin >> dim;
  Linear::Matrix A(dim, dim);
  Linear::Vector v(dim), res(dim);

  A.Fill(std::cin);
  v.Fill(std::cin);

  start = omp_get_wtime();
  Linear::Mult(A, v, res);
  end = omp_get_wtime();
  
  std::cout << end - start << std::endl;
}
