/*
 * The parallel calculation of the product of two square matrices.
 */
#include <iostream>
#include "linear.h"
#include <omp.h>

int main(void) {
  int thNum, dim;
  double start, end;

  std::cin >> thNum >> dim;
  Linear::Matrix A(dim, dim), B(dim, dim), res(dim, dim);
  A.Fill(std::cin);
  B.Fill(std::cin);

  start = omp_get_wtime();
  Linear::PMult(A, B, res, thNum);
  end = omp_get_wtime();
  
  std::cout << end - start << std::endl;
}
