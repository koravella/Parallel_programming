/*
 * Calculates the LU decomposition of square matrix.
 */
#include <iostream>
#include "linear.h"
#include <omp.h>

int main(void) {
  int dim;
  double start, end;

  std::cin >> dim;
  Linear::Matrix A(dim, dim), L(dim, dim), U(dim, dim);
  A.Fill(std::cin);

  start = omp_get_wtime();
  Linear::LUdecomposition(A, L, U);
  end = omp_get_wtime();
  
  std::cout << end - start << std::endl;
  //L.Print(std::cout);
  //U.Print(std::cout);
}
