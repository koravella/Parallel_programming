#include "linear.h"
#include <iostream>
#include <cmath>
#include <omp.h>

using namespace Linear;

Vector::Vector(const int n) {
  this->size = n;
  this->data = new double[n];
}

Vector::~Vector() {
  delete[] this->data;
}

double& Vector::operator()(const int index) {
  return data[index];
}

double Vector::operator()(const int index) const {
  return data[index];
}

int Vector::Getsize() const {
  return this->size;
}

void Vector::Fill(std::istream& fin) {
  for (int i = 0; i < this->size; i++)
    fin >> this->data[i];
}

void Vector::Print(std::ostream& fout) const {
  for (int i = 0; i < this->size; i++)
      fout << this->data[i] << " ";
  fout << std::endl;
}

Matrix::Matrix(const int row, const int col) {
  this->nrow = row;
  this->ncol = col;
  this->data = new double[row * col];
}

Matrix::~Matrix() {
  delete[] this->data;
}

double& Matrix::operator()(const int indrow, const int indcol) {
  return data[indrow * ncol + indcol];
}

double Matrix::operator()(const int indrow, const int indcol) const {
  return data[indrow * ncol + indcol];
}

int Matrix::Getnrow() const {
  return this->nrow;
}

int Matrix::Getncol() const {
  return this->ncol;
}

void Matrix::Fill(std::istream& fin) {
  for (int i = 0; i < this->nrow; i++)
    for (int k = 0; k < this->ncol; k++)
      fin >> this->data[i * this->ncol + k];
}

void Matrix::Print(std::ostream& fout) const {
  for (int i = 0; i < this->nrow; i++) {
    for (int k = 0; k < this->ncol; k++)
      fout << this->data[i * this->ncol + k] << " ";
    fout << std::endl;
  }
}

int Linear::Mult(const Matrix& A, const Matrix& B, Matrix& res) {
  if (A.Getncol() != B.Getnrow() || res.Getnrow() != A.Getnrow()
      || res.Getncol() != B.Getncol())
    return -1; // error
  
  for (int row = 0; row < res.Getnrow(); row++)
    for (int col = 0; col < res.Getncol(); col++) {
      res(row, col) = 0;
      for (int k = 0; k < A.Getncol(); k++)
        res(row, col) += A(row, k) * B(k, col);
    }
  return 0;
}

int Linear::Mult(const Matrix& A, const Vector& v, Vector& res) {
  if (A.Getncol() != v.Getsize() || res.Getsize() != A.Getnrow())
    return -1; // error

  for (int cur = 0; cur < res.Getsize(); cur++) {
    res(cur) = 0;
    for (int k = 0; k < A.Getncol(); k++)
      res(cur) += A(cur, k) * v(k);
  }
  return 0;
}

int Linear::LUdecomposition(const Matrix& A, Matrix& L, Matrix& U) {
  int n;

  if ((n = A.Getnrow()) != A.Getncol() ||
      n != L.Getnrow() || n != L.Getncol() ||
      n != U.Getnrow() || n != U.Getncol())
    return -1; // error
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (j < i) {
        L(j, i) = 0;
        continue;
      }
      L(j, i) = A(j, i);
      for (int k = 0; k < i; k++)
        L(j, i) -= L(j, k) * U(k, i);
    }
    for (int j = 0; j < n; j++) {
      if (j < i) {
        U(i, j) = 0;
        continue;
      }
      if (j == i) {
        U(i, j) = 1;
        continue;
      }
//      if (abs(L(i,i)) < eps)
//        return -1; // LU does not exist
      U(i, j) = A(i, j) / L(i, i);
      for (int k = 0; k < i; k++)
        U(i, j) -= L(i, k) * U(k, j) / L(i, i);
      //U(i, j) /= L(i, i);
    }
  }
  return 0;
}

int Linear::Inversion(const Matrix& A, Matrix& res) {
  int n;
  double det = 1;

  if ((n = A.Getnrow()) != A.Getncol() ||
      n != res.Getnrow() || n != res.Getncol())
    return -1; // error
  
  Matrix L(n, n), U(n, n);
  int ret = LUdecomposition(A, L, U);
//  if (ret == -1)
//    return -1; // LU does not exist
  
  for (int i = 0; i < n; i++) {
    if (abs(L(i, i)) < eps || abs(U(i, i)) < eps) {
      det = 0;
      break;
    }
    det *= L(i, i) * U(i, i);
  }
  if (abs(det) < eps)
    return -2;  // determinant equals zero

  for (int t = 0; t < n; t++) { // columns of res
    Vector y(n);
    
    for (int k = 0; k < n; k++) { // step 1 forward
      double sum = k == t ? 1 : 0; // identity matrix column
      for (int i = 0; i < k; i++)
        sum -= L(k, i) * y(i);
      y(k) = sum / L(k, k);
    }
    
    for (int k = n - 1; k >= 0; k--) { // step 2 backward
      double sum = y(k);
      for (int i = n - 1; i > k; i--)
        sum -= U(k, i) * res(i, t);
      res(k, t) = sum / U(k, k);
    }
  }
  return 0;
}

int Linear::PMult(const Matrix& A, const Matrix& B, Matrix& res, int thNum) {
  if (A.Getncol() != B.Getnrow() || res.Getnrow() != A.Getnrow()
      || res.Getncol() != B.Getncol())
    return -1; // error

  int row, col, k;
  omp_set_num_threads(thNum);
#pragma omp parallel for shared(A, B, res) private(row, col, k)
  for (row = 0; row < res.Getnrow(); row++)
    for (col = 0; col < res.Getncol(); col++) {
      res(row, col) = 0;
      for (k = 0; k < A.Getncol(); k++)
        res(row, col) += A(row, k) * B(k, col);
    }
  return 0;
}

int Linear::PMult(const Matrix& A, const Vector& v, Vector& res, int thNum) {
  if (A.Getncol() != v.Getsize() || res.Getsize() != A.Getnrow())
    return -1; // error

  int cur, k;
  omp_set_num_threads(thNum);
#pragma omp parallel for shared(A, v, res) private(cur, k) 
  for (cur = 0; cur < res.Getsize(); cur++) {
    res(cur) = 0;
    for (k = 0; k < A.Getncol(); k++)
      res(cur) += A(cur, k) * v(k);
  }
  return 0;
}

int Linear::PLUdecomposition(const Matrix& A, Matrix& L, Matrix& U, int thNum) {
  int n;

  if ((n = A.Getnrow()) != A.Getncol() ||
      n != L.Getnrow() || n != L.Getncol() ||
      n != U.Getnrow() || n != U.Getncol())
    return -1; // error

  int flag = 0;
#pragma omp parallel shared(A, L, U, flag)
  {
//  #pragma omp for
  for (int i = 0; i < n; i++) {

#pragma omp for schedule(dynamic, 32)
    for (int j = 0; j < n; j++) {
      if (j < i) {
        L(j, i) = 0;
        continue;
      }
      L(j, i) = A(j, i);
      for (int k = 0; k < i; k++)
        L(j, i) -= L(j, k) * U(k, i);
    }
#pragma omp for schedule(dynamic, 32)
    for (int j = 0; j < n; j++) {
      if (j < i) {
        U(i, j) = 0;
        continue;
      }
      if (j == i) {
        U(i, j) = 1;
        continue;
      }
      if (abs(L(i, i) < eps)) {
#pragma omp critical
        {
        flag = 1; // LU does not exist
        }
      }
      U(i, j) = A(i, j);
      for (int k = 0; k < i; k++)
        U(i, j) -= L(i, k) * U(k, j);
      U(i, j) /= L(i, i);
    }
  }
  }
//  if (flag == 1)
//    return -1;
  return 0;
}

int Linear::PInversion(const Matrix& A, Matrix& res, int thNum) {
  int n;
  double det = 1;

  if ((n = A.Getnrow()) != A.Getncol() ||
      n != res.Getnrow() || n != res.Getncol())
    return -1; // error
  
  Matrix L(n, n), U(n, n);
  int ret = PLUdecomposition(A, L, U, thNum);
//  if (ret == -1)
//    return -1; // LU does not exist
  
  int i;
#pragma omp parallel for shared(L, U) private(i) reduction(*:det)
  for (i = 0; i < n; i++) {
    if (abs(L(i, i)) < eps || abs(U(i, i)) < eps)
      det *= 0;
    det *= L(i, i) * U(i, i);
  }
  if (abs(det) < eps)
    return -2;  // determinant equals zero

#pragma omp parallel shared(L, U, res)
  {
#pragma omp for
  for (int t = 0; t < n; t++) { // columns of res
    Vector y(n);
    
    for (int k = 0; k < n; k++) { // step 1 forward
      double sum = k == t ? 1 : 0; // identity matrix column
      for (int i = 0; i < k; i++)
        sum -= L(k, i) * y(i);
      y(k) = sum / L(k, k);
    }
    
    for (int k = n - 1; k >= 0; k--) { // step 2 backward
      double sum = y(k);
      for (int i = n - 1; i > k; i--)
        sum -= U(k, i) * res(i, t);
      res(k, t) = sum / U(k, k);
    }
  }
  }
  return 0;
}

