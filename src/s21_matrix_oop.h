#ifndef CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
#define CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H

#include <cmath>
#include <exception>
#include <iostream>

namespace Matrix {

class S21Matrix {
 public:
  // constructors and destructor
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();

  // methods
  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(double num);
  void MulMatrix(const S21Matrix& other);
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix Transpose();
  S21Matrix InverseMatrix();
  int GetRows() const;
  int GetCols() const;
  void SetRows(int count);
  void SetCols(int count);

  // operators
  S21Matrix& operator=(const S21Matrix& other) noexcept;  // copy
  S21Matrix& operator=(S21Matrix&& other) noexcept;       // move
  bool operator==(const S21Matrix& other) noexcept;
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(double num);
  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  void operator*=(double num);
  double& operator()(int row, int col);
  double operator()(int row, int col) const;

 private:
  int rows_ = 0, cols_ = 0;
  double** matrix_ = nullptr;

  // my addings
  void MemoryAlloc();
  void SumSubMatrix(const S21Matrix& other, int sign);
  S21Matrix TakeMinor(int rows, int cols);
  void RemoveMatrix();
  void MatrixCheck() const;
  void SwapMatrix(S21Matrix& other) noexcept;
};

}  // namespace S21MatrixNames

#endif  //  CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
