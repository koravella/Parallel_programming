/*
 * The parallel calculation of the inversion of a square matrix.
 */
#include <iostream>
#include "linear.h"
#include <omp.h>

int main(void) {
  int thNum, dim;
  double start, end;

  std::cin >> thNum >> dim;
  Linear::Matrix A(dim, dim), res(dim, dim);
  A.Fill(std::cin);

  start = omp_get_wtime();
  Linear::PInversion(A, res, thNum);
  end = omp_get_wtime();
  
  std::cout << end - start << std::endl;
}
