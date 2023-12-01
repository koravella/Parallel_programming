#ifndef LINEAR_H
#define LINEAR_H

#include <fstream>

#define eps 0.000001  // small value 

namespace Linear {
  
  class Vector {
    int size;
    double *data;
  public:
    Vector() = delete;
    Vector(const int n);
    ~Vector();
    Vector(const Vector& that) = delete;
    Vector(Vector&& that) = delete;
    Vector& operator=(const Vector& that) = delete;
    double& operator()(const int index);
    double operator()(const int index) const;
    int Getsize() const;
    void Fill(std::istream& fin);
    void Print(std::ostream& fout) const;
  };
  
  class Matrix {
    int nrow;
    int ncol;
    double *data;
  public:
    Matrix() = delete;
    Matrix(const int row, const int col);
    ~Matrix();
    Matrix(const Matrix& that) = delete;
    Matrix(Matrix&& that) = delete;
    Matrix& operator=(const Matrix& that) = delete;
    double& operator()(const int indrow, const int indcol);
    double operator()(const int indrow, const int indcol) const;
    int Getnrow() const;
    int Getncol() const;
    void Fill(std::istream& fin);
    void Print(std::ostream& fout) const;
  };
  
  /* 
   * The product of two matrices A and B: C = A * B.
   * Return -1 if error and 0 otherwise.
   */
  int Mult(const Matrix& A, const Matrix& B, Matrix& res);
  
  /*
   * The product of a matrix A and vector v: res = A * v.
   * Return -1 if error and 0 otherwise.
   */
  int Mult(const Matrix& A, const Vector& v, Vector& res);
  
  /*
   * LU decomposition of matrix A: A = L * U, where L is lower triangular and
   * U is upper triangular.
   * Return -1 if error and 0 otherwise.
   */
  int LUdecomposition(const Matrix& A, Matrix& L, Matrix& U);
  
  /*
   * Inversion of matrix A: res = A^(-1).
   * Return -2 if determinant is equal zero, -1 if error and 0 otherwise.
   */
  int Inversion(const Matrix& A, Matrix& res);

  /* 
   * The parallel calculation of product of two matrices A and B: C = A * B.
   * thNum -- number of threads.
   * Return -1 if error and 0 otherwise.
   */
  int PMult(const Matrix& A, const Matrix& B, Matrix& res, int thNum);
  
  /*
   * The parallel calculation of the product of a matrix A and vector v: res = A * v.
   * thNum -- number of threads.
   * Return -1 if error and 0 otherwise.
   */
  int PMult(const Matrix& A, const Vector& v, Vector& res, int thNum);

  /*
   * Parallel LU decomposition of matrix A: A = L * U, where L is lower
   * triangular and U is upper triangular.
   * thNum -- number of threads.
   * Return -1 if error and 0 otherwise.
   */
  int PLUdecomposition(const Matrix& A, Matrix& L, Matrix& U, int thNum);
   
  /*
   * Parallel inversion of matrix A: res = A^(-1).
   * thNum -- number of threads.
   * Return -2 if determinant is equal zero, -1 if error and 0 otherwise.
   */
  int PInversion(const Matrix& A, Matrix& res, int thNum);
}



#endif
