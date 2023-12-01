/*
 * Calculates the product of two square matrices.
 */
#include <iostream>
#include "linear.h"
#include <omp.h>

int main(void) {
  int dim;
  double start, end;

  std::cin >> dim;
  Linear::Matrix A(dim, dim), B(dim, dim), res(dim, dim);
  A.Fill(std::cin);
  B.Fill(std::cin);

  start = omp_get_wtime();
  Linear::Mult(A, B, res);
  end = omp_get_wtime();
  
  std::cout << end - start << std::endl;
}
