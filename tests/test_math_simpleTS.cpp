#include "math.h"
#include "matrix.h"

int main(int argc, char* argv[]) {
  CHECK(argc == 3, "Usage: ./a.out A_file.mtx B_file.mtx");

  std::string A_file(argv[1]), b_file(argv[2]);
  std::unique_ptr<CSSparseMatrix<double>> A =
      std::move(CSSparseMatrix<double>::ConstructFromFile(A_file));
  std::unique_ptr<DenseMatrix<double>> b = std::move(
      DenseMatrix<double>::ConstructFromFile(b_file, DenseStorageType::RMaj));

  {
    std::unique_ptr<DenseMatrix<double>> b_clone = std::move(b->clone());
    math::SimpleSparseTriangularSolve<double>(A.get(), b_clone.get());
    for (size_t i = 0; i < b_clone->rows(); ++i) {
      for (size_t j = 0; j < b_clone->cols(); ++j) {
        std::cout << *(b_clone->access(i, j)) << ' ';
      }
      std::cout << '\n';
    }
  }

  {
    std::unique_ptr<DenseMatrix<double>> b_clone = std::move(b->clone());
    math::Opt1SparseTriangularSolve<double>(A.get(), b_clone.get());
    for (size_t i = 0; i < b_clone->rows(); ++i) {
      for (size_t j = 0; j < b_clone->cols(); ++j) {
        std::cout << *(b_clone->access(i, j)) << ' ';
      }
      std::cout << '\n';
    }
  }
  return 0;
}
