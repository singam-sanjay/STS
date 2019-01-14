#include "math.h"
#include "matrix.h"

int main(int argc, char* argv[]) {
  CHECK(argc == 3, "Usage: ./a.out A_file.mtx B_file.mtx");

  std::string A_file(argv[1]), b_file(argv[2]);
  std::unique_ptr<CSSparseMatrix<double>> A =
      std::move(CSSparseMatrix<double>::ConstructFromFile(A_file));
  std::unique_ptr<DenseMatrix<double>> b = std::move(
      DenseMatrix<double>::ConstructFromFile(b_file, DenseStorageType::RMaj));

  math::SimpleSparseTriangularSolve<double>(A.get(), b.get());
  for (size_t i = 0; i < b->rows(); ++i) {
    for (size_t j = 0; j < b->cols(); ++j) {
      std::cout << *(b->access(i, j)) << ' ';
    }
    std::cout << '\n';
  }
  return 0;
}
