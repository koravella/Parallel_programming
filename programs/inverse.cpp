/*
 * Calculates the inversion of a square matrix.
 */
#include <iostream>
#include "linear.h"
#include <omp.h>

int main(void) {
  int dim;
  double start, end;

  std::cin >> dim;
  Linear::Matrix A(dim, dim), res(dim, dim);
  A.Fill(std::cin);

  start = omp_get_wtime();
  Linear::Inversion(A, res);
  end = omp_get_wtime();
  
  std::cout << end - start << std::endl;
  //res.Print(std::cout);
}
