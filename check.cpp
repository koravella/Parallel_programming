#include "linear.h"
#include <iostream>
#include <omp.h>

#define thNum 8

int main(void) {
  int r1, c1, r2, c2, r3;

  double start = omp_get_wtime();

  std::cin >> r1 >> c1;
  Linear::Matrix A(r1, c1);
  A.Fill(std::cin);

  std::cin >> r2 >> c2;
  Linear::Matrix B(r2, c2);
  B.Fill(std::cin);
  
  Linear::Matrix res(r1, c2), res_par(r1, c2);
  Linear::Mult(A, B, res);
  std::cout <<"A * B ="<< std::endl;
  res.Print(std::cout);
  Linear::PMult(A, B, res_par, thNum);
  std::cout <<"\n[parallel]\nA * B ="<< std::endl;
  res_par.Print(std::cout);

  std::cout << "\nA * v =" << std::endl;
  Linear::Vector v(c1), resv(r1), resv_par(r1);
  v.Fill(std::cin);
  Linear::Mult(A, v, resv);
  resv.Print(std::cout);
  Linear::PMult(A, v, resv_par, thNum);
  std::cout << "\n[parallel]\nA * v =" << std::endl;
  resv_par.Print(std::cout);

  std::cin >> r3;
  Linear::Matrix C(r3,r3), L(r3,r3), U(r3,r3), D(r3,r3), L_par(r3,r3), U_par(r3,r3);
  C.Fill(std::cin);
  Linear::LUdecomposition(C, L, U);
  std::cout << "\nL =" << std::endl;
  L.Print(std::cout);
  std::cout << "\nU =" << std::endl;
  U.Print(std::cout);
  Linear::Mult(L, U, D);
  std::cout << "\nD = L * U =" << std::endl;
  D.Print(std::cout);
  Linear::PLUdecomposition(C, L_par, U_par, thNum);
  std::cout << "\n[parallel]\nL =" << std::endl;
  L_par.Print(std::cout);
  std::cout << "\nU =" << std::endl;
  U_par.Print(std::cout);


  std::cout << "\nC^-1 = " << std::endl;
  Linear::Matrix E(r3, r3), E_par(r3,r3);
  if (Linear::Inversion(C, E) == -2)
    std::cout << "Determinant is equal zero\n";
  else
    E.Print(std::cout);
  if (Linear::PInversion(C, E_par, thNum) == -2)
    std::cout << "\n[parallel]\nDeterminant is equal zero\n";
  else {
    std::cout <<"\n[parallel]\n";
    E_par.Print(std::cout);
  }


  double end = omp_get_wtime();
  std::cout << "\nTime: " << end - start << std::endl;

}
