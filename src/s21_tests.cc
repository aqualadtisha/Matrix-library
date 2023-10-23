#include <gtest/gtest.h>

#include <vector>

#include "s21_matrix_oop.h"

namespace S21MatrixNames {

class TestS21Matrix : public ::testing::Test {
 public:
  //  void DataFill(S21Matrix m, const double* Data);

 protected:
  S21Matrix* matrix1;
  S21Matrix* matrix3;
  S21Matrix* matrix5;
  S21Matrix* matrix7;

  double data[4] = {2, 1, 2, 4};
  double dataout[4] = {2, 0, 2, 4};
  double dataout2[4] = {4, 2, 4, 8};
  double dataout4[4] = {6, 6, 12, 18};

  //  void PrintMatrix(S21Matrix m);
  //  TestS21Matrix() {}

  void SetUp() {
    matrix1 = new S21Matrix(2, 2);
    DataFill(*matrix1, data);
    matrix3 = new S21Matrix(2, 2);
    DataFill(*matrix3, dataout);
    matrix5 = new S21Matrix(2, 2);
    DataFill(*matrix5, dataout2);
    matrix7 = new S21Matrix(2, 2);
    DataFill(*matrix7, dataout4);
  }
  void TearDown() {
    delete matrix1;
    delete matrix3;
    delete matrix5;
    delete matrix7;
  }

  void DataFill(S21Matrix& m, const double* Data) {
    int k = 0;
    int rows = m.GetRows();
    int cols = m.GetCols();
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++) m(i, j) = Data[k++];
  }
  void PrintMatrix(S21Matrix m) {
    int rows = m.GetRows();
    int cols = m.GetCols();
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) printf("%f ", m(i, j));
      printf("\n");
    }
    printf("\n");
  }
};

TEST_F(TestS21Matrix, CreateMatrix) {
  S21Matrix mat(2, 2);
  EXPECT_EQ(mat.GetCols(), 2);
  EXPECT_THROW(S21Matrix n(-1, 2), std::exception);
  EXPECT_THROW(S21Matrix n(1, 0), std::exception);
}

TEST_F(TestS21Matrix, MoveMatrix) {
  S21Matrix matrix2(std::move(*matrix1));
  EXPECT_EQ(matrix2.GetRows() == 2, true);
  EXPECT_EQ(matrix2.GetCols() == 2, true);
}

TEST_F(TestS21Matrix, GetRowsColsMatrix) {
  S21Matrix mat(2, 2);
  int rows = mat.GetRows();
  EXPECT_EQ(rows, 2);
  EXPECT_EQ(mat.GetCols(), 2);
  EXPECT_THROW(mat.SetRows(0), std::exception);
}

TEST_F(TestS21Matrix, SetRowsColsMatrix) {
  S21Matrix mat(2, 2);
  mat.SetCols(5);
  mat.SetRows(6);
  int rows = mat.GetRows();
  EXPECT_EQ(rows, 6);
  EXPECT_EQ(mat.GetCols(), 5);
  EXPECT_THROW(mat.SetCols(0), std::exception);
}
TEST_F(TestS21Matrix, EqMatrix) {
  S21Matrix matrix2(*matrix1);

  EXPECT_EQ(matrix1->EqMatrix(matrix2), true);
  EXPECT_EQ(matrix1->EqMatrix(*matrix3), false);
}

TEST_F(TestS21Matrix, SumMatrixFunction) {
  S21Matrix matrix2(*matrix1);

  matrix1->SumMatrix(matrix2);

  EXPECT_EQ(matrix1->EqMatrix(*matrix5), true);
  matrix1->SetCols(7);
  EXPECT_THROW(*matrix1 + *matrix5, std::exception);
}

TEST_F(TestS21Matrix, SubMatrixFunction) {
  S21Matrix matrix2(*matrix1);
  S21Matrix matrix3(2, 2);

  matrix1->SubMatrix(matrix2);
  EXPECT_EQ(matrix1->EqMatrix(matrix3), true);
  matrix1->SetCols(7);
  EXPECT_THROW(*matrix1 - matrix3, std::exception);
}

TEST_F(TestS21Matrix, MulNumberFunction) {
  double dataout[4] = {8, 4, 8, 16};

  S21Matrix matrix3(2, 2);
  DataFill(matrix3, dataout);
  matrix1->MulNumber(4);
  EXPECT_EQ(matrix1->EqMatrix(matrix3), true);
}
TEST_F(TestS21Matrix, MulMatrixFunction) {
  S21Matrix matrix2(*matrix1);
  matrix1->MulMatrix(matrix2);
  EXPECT_EQ(matrix1->EqMatrix(*matrix7), true);
}

TEST_F(TestS21Matrix, MulMatrixFunctionWithExeption) {
  double dataout[6] = {2, 1, 2, 4, 0, 0};

  S21Matrix matrix3(3, 2);
  DataFill(matrix3, dataout);
  matrix1->SetRows(3);
  EXPECT_EQ(matrix1->EqMatrix(matrix3), true);
  EXPECT_THROW(*matrix1 * matrix3, std::exception);
}

TEST_F(TestS21Matrix, TrancFunction) {
  double dataout[4] = {2, 2, 1, 4};

  S21Matrix matrix2(*matrix1);
  matrix2 = matrix1->Transpose();
  S21Matrix matrix3(2, 2);
  DataFill(matrix3, dataout);
  EXPECT_EQ(matrix2.EqMatrix(matrix3), true);
}
TEST_F(TestS21Matrix, CalcComplemencFunction) {
  double dataout[4] = {4, -2, -1, 2};

  S21Matrix matrix2;
  matrix2 = matrix1->CalcComplements();
  S21Matrix matrix3(2, 2);
  DataFill(matrix3, dataout);
  EXPECT_EQ(matrix2.EqMatrix(matrix3), true);
  matrix2.SetRows(4);
  EXPECT_THROW(matrix2.CalcComplements(), std::exception);
  S21Matrix matrixNULL = std::move(matrix2);
  EXPECT_THROW(matrix2.CalcComplements(), std::exception);

  S21Matrix matrix4(1, 1);
  double d[1] = {2};
  DataFill(matrix4, d);
  EXPECT_THROW(matrix4.CalcComplements(), std::exception);
}

TEST_F(TestS21Matrix, DeterminantFunction) {
  double det = matrix1->Determinant();
  EXPECT_EQ(det, 6);
  matrix1->SetRows(4);
  EXPECT_THROW(matrix1->Determinant(), std::exception);

  double data1[9] = {1, 2, 33, 2, 543, 23, 123, 123, 124};
  S21Matrix matrix2(3, 3);
  DataFill(matrix2, data1);
  EXPECT_EQ(matrix2.Determinant() == -2126254, true);

  double data2[16] = {1, 2, 3, 4, 5, -6, 7, 8, 9, 10, -11, -12, 13, 14, 15, 16};
  S21Matrix matrix(4, 4);
  DataFill(matrix, data2);
  EXPECT_EQ(matrix.Determinant() == -2592, true);
}

TEST_F(TestS21Matrix, InverseMatrixFunction) {
  double dataout[4] = {2 / 3., -1 / 6., -1 / 3., 1 / 3.};
  double data1[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  S21Matrix matrix2(*matrix1);
  matrix2 = matrix1->InverseMatrix();
  S21Matrix matrix3(2, 2);
  DataFill(matrix3, dataout);
  EXPECT_EQ(matrix2.EqMatrix(matrix3), true);
  matrix3.SetRows(8);
  EXPECT_THROW(matrix3.InverseMatrix(), std::exception);
  S21Matrix ZeroDet(3, 3);
  DataFill(ZeroDet, data1);
  EXPECT_THROW(ZeroDet.InverseMatrix(), std::exception);

  double data2[1] = {2};

  S21Matrix matrix4(1, 1);
  DataFill(matrix4, data2);
  S21Matrix tt = matrix4.InverseMatrix();
  EXPECT_EQ(tt(0, 0) == 0.5, true);
}

TEST_F(TestS21Matrix, EqualToOperator) {
  S21Matrix matrix2(*matrix1);

  *matrix1 = matrix2 + *matrix1;
  EXPECT_EQ((*matrix1 == *matrix5), true);
  EXPECT_EQ((*matrix1 == matrix2), false);
}

TEST_F(TestS21Matrix, SumMatrixOperator) {
  S21Matrix matrix2(*matrix1);
  *matrix1 = matrix2 + *matrix1;
  EXPECT_EQ(matrix1->EqMatrix(*matrix5), true);
}

TEST_F(TestS21Matrix, SumMatrixAssignmentOperator) {
  S21Matrix matrix2(*matrix1);
  *matrix1 += matrix2;
  EXPECT_EQ(matrix1->EqMatrix(*matrix5), true);
}

TEST_F(TestS21Matrix, SubMatrixOperator) {
  S21Matrix matrix2(*matrix1);
  S21Matrix matrix4(2, 2);
  *matrix1 = matrix2 - *matrix1;
  EXPECT_EQ(matrix1->EqMatrix(matrix4), true);
}

TEST_F(TestS21Matrix, SubMatrixAssignmentOperator) {
  S21Matrix matrix2(*matrix1);
  S21Matrix matrix3(2, 2);
  *matrix1 -= matrix2;
  EXPECT_EQ(matrix1->EqMatrix(matrix3), true);
}

TEST_F(TestS21Matrix, MulMatrixOperator) {
  S21Matrix matrix2(*matrix1);

  *matrix1 = matrix2 * *matrix1;
  EXPECT_EQ(matrix1->EqMatrix(*matrix7), true);
}

TEST_F(TestS21Matrix, MulMatrixAssignmentOperator) {
  S21Matrix matrix2(*matrix1);

  *matrix1 *= matrix2;
  EXPECT_EQ(matrix1->EqMatrix(*matrix7), true);
}

TEST_F(TestS21Matrix, MulNumberOperator) {
  *matrix1 = *matrix1 * 2;
  EXPECT_EQ(matrix5->EqMatrix(*matrix1), true);
}
TEST_F(TestS21Matrix, MulNumberFriendOperator) {
  *matrix1 = *matrix1 * 2.0;
  EXPECT_EQ(matrix5->EqMatrix(*matrix1), true);
}

TEST_F(TestS21Matrix, MulNumberAssignmentOperator) {
  *matrix1 *= 2;
  EXPECT_EQ(matrix5->EqMatrix(*matrix1), true);
}
TEST_F(TestS21Matrix, OpertorGetElement) {
  EXPECT_EQ(matrix1->operator()(0, 0) == data[0], true);
  EXPECT_THROW(matrix1->operator()(-1, 2), std::exception);
  EXPECT_THROW(matrix1->operator()(3, 2), std::exception);
}

}  // namespace S21MatrixNames

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}