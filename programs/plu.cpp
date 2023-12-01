/*
 * The parallel calculation of the LU decomposition of square matrix.
 */
#include <iostream>
#include "linear.h"
#include <omp.h>

int main(void) {
  int thNum, dim;
  double start, end;

  std::cin >> thNum >> dim;
  Linear::Matrix A(dim, dim), L(dim, dim), U(dim, dim);
  A.Fill(std::cin);

  start = omp_get_wtime();
  Linear::PLUdecomposition(A, L, U, thNum);
  end = omp_get_wtime();
  
  std::cout << end - start << std::endl;
}
