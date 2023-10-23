#include "s21_matrix_oop.h"

namespace Matrix {

S21Matrix::S21Matrix() = default;

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) throw std::range_error("Invalid matrix size");
  cols_ = cols;
  rows_ = rows;
  MemoryAlloc();
}

void S21Matrix::MemoryAlloc() {
  matrix_ = new double*[rows_]();
  for (auto i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

S21Matrix::S21Matrix(const S21Matrix& other) {  // copy
  rows_ = other.rows_;
  cols_ = other.cols_;
  MemoryAlloc();

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other(i, j);
    }
  }
}


S21Matrix::S21Matrix(S21Matrix&& other) noexcept {  // move
  cols_ = std::move(other.cols_);
  rows_ = std::move(other.rows_);
  matrix_ = std::move(other.matrix_);
}

S21Matrix::~S21Matrix() { RemoveMatrix(); }

void S21Matrix::RemoveMatrix() {
  if (matrix_) {
    for (auto i = 0; i < rows_; i++) delete[] matrix_[i];
    delete[] matrix_;
  }

  rows_ = 0;
  cols_ = 0;
  matrix_ = nullptr;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  if (this == &other) return true;

  if (rows_ != other.rows_ || cols_ != other.cols_) return false;
  if (matrix_ == nullptr && other.matrix_ == nullptr) return true;

  bool error = true;

  for (int i = 0; i < rows_ && error; i++) {
    for (int j = 0; j < cols_ && error; j++) {
      if (matrix_[i] && other.matrix_[i]) {
        if (fabs(matrix_[i][j] - other(i, j)) >= 1e-7) error = false;
      } else {
        error = false;
      }
    }
  }

  return error;
}

void S21Matrix::SumMatrix(const S21Matrix& other) { SumSubMatrix(other, 1); }
void S21Matrix::SubMatrix(const S21Matrix& other) { SumSubMatrix(other, -1); }

void S21Matrix::SumSubMatrix(const S21Matrix& other, int sign) {
  MatrixCheck();
  other.MatrixCheck();

  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::length_error("Invalid second matrix size");

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] + (sign * other(i, j));
    }
  }
}

void S21Matrix::MatrixCheck() const {
  if (matrix_ == nullptr || rows_ <= 0 || cols_ <= 0)
    throw std::invalid_argument("Invalid matrix");
}

void S21Matrix::MulNumber(const double num) {
  MatrixCheck();

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  MatrixCheck();
  other.MatrixCheck();

  if (cols_ != other.rows_)
    throw std::length_error("Invalid second matrix size");

  S21Matrix tmp(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int k = 0; k < other.cols_; k++) {
      for (int side = 0; side < cols_; side++) {
        tmp.matrix_[i][k] += matrix_[i][side] * other(side, k);
      }
    }
  }

  *this = tmp;
}

double S21Matrix::Determinant() {
  double result = 0;
  MatrixCheck();

  if (rows_ != cols_) throw std::range_error("Matrix isn't square");

  if (rows_ < 3) {

    if (rows_ == 2) {
      result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
    } else {
      result = matrix_[0][0];
    }

  } else {

    S21Matrix tmp(rows_ - 1, cols_ - 1);

    for (int column = 0; column < cols_; column++) {
      int sign = column % 2 == 0 ? 1 : -1;
      tmp = TakeMinor(0, column);

      if (tmp.rows_ == 2) {
        result += sign * matrix_[0][column] *
                  (tmp.matrix_[0][0] * tmp.matrix_[1][1] -
                   tmp.matrix_[0][1] * tmp.matrix_[1][0]);
      } else {
        double add_res = tmp.Determinant();
        result += sign * matrix_[0][column] * add_res;
      }

    }

  }
  return result;
}

S21Matrix S21Matrix::TakeMinor(int rows, int cols) {
  S21Matrix tmp(cols_ - 1, rows_ - 1);
  int k = 0;

  for (int i = 0; i < rows_; i++) {
    int m = 0;
    if (i != rows) {
      for (int j = 0; j < cols_; j++) {
        if (j != cols) {
          tmp.matrix_[k][m++] = matrix_[i][j];
        }
      }
      k++;
    }
  }

  return tmp;
}

S21Matrix S21Matrix::CalcComplements() {
  MatrixCheck();
  if (rows_ != cols_) throw std::range_error("Matrix isn't square");

  S21Matrix result(rows_, cols_);
  S21Matrix tmp;

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < rows_; j++) {

      int sign = (j + i) % 2 == 0 ? 1 : -1;
      tmp = TakeMinor(i, j);
      result.matrix_[i][j] = sign * tmp.Determinant();
    }
  }

  return result;
}

S21Matrix S21Matrix::Transpose() {
  MatrixCheck();
  S21Matrix tmp(cols_, rows_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      tmp.matrix_[j][i] = matrix_[i][j];
    }
  }

  return tmp;
}

S21Matrix S21Matrix::InverseMatrix() {
  MatrixCheck();
  if (rows_ != cols_) throw std::range_error("Matrix isn't square");

  double det = Determinant();
  if (det == 0) throw std::invalid_argument("Zero determinant");

  S21Matrix result(1, 1);

  if (rows_ == 1) {
    result.matrix_[0][0] = 1 / det;
  } else {
    result = CalcComplements().Transpose() * (1 / det);
  }

  return result;
}

S21Matrix& S21Matrix::operator=(
    S21Matrix&& other) noexcept {  // move assignment
  if (this == &other) return *this;

  S21Matrix tmp((std::move(other)));
  SwapMatrix(tmp);

  return *this;
}

S21Matrix& S21Matrix::operator=(
    const S21Matrix& other) noexcept {  // copy assignment
  if (this == &other) return *this;

  S21Matrix tmp(other);
  SwapMatrix(tmp);

  return *this;
}

void S21Matrix::SwapMatrix(S21Matrix& other) noexcept {
  std::swap(cols_, other.cols_);
  std::swap(rows_, other.rows_);
  std::swap(matrix_, other.matrix_);
}

bool S21Matrix::operator==(const S21Matrix& other) noexcept {
  return EqMatrix(other);
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(double num) {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}

void S21Matrix::operator+=(const S21Matrix& other) { SumMatrix(other); }

void S21Matrix::operator-=(const S21Matrix& other) { SubMatrix(other); }

void S21Matrix::operator*=(const S21Matrix& other) { MulMatrix(other); }

void S21Matrix::operator*=(double num) { MulNumber(num); }

double& S21Matrix::operator()(int row, int col) {
  if (row >= rows_ || row < 0 || col < 0 || col >= cols_)
    throw std::out_of_range("Index out of range");

  return matrix_[row][col];
}

double S21Matrix::operator()(int row, int col) const {
  if (row >= rows_ || row < 0 || col < 0 || col >= cols_)
    throw std::out_of_range("Index out of range");

  return matrix_[row][col];
}

int S21Matrix::GetRows() const { return rows_; }
int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetRows(int count) {
  S21Matrix tmp(count, cols_);

  for (int i = 0; i < count && i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }

  *this = tmp;
}

void S21Matrix::SetCols(int count) {
  S21Matrix tmp(rows_, count);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < count && j < cols_; j++) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }

  *this = tmp;
}

}  // namespace S21MatrixNames