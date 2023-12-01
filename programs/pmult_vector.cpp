/*
 * Parallel calculation of the product of square matrix and vector.
 */
#include <iostream>
#include "linear.h"
#include <omp.h>

int main(void) {
  int thNum, dim;
  double start, end;

  std::cin >> thNum >> dim;
  Linear::Matrix A(dim, dim);
  Linear::Vector v(dim), res(dim);

  A.Fill(std::cin);
  v.Fill(std::cin);

  start = omp_get_wtime();
  Linear::PMult(A, v, res, thNum);
  end = omp_get_wtime();
  
  std::cout << end - start << std::endl;
}
