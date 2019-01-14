#include <chrono>

#include "math.h"
#include "matrix.h"

#define TICK() auto ____NOW____ = std::chrono::system_clock::now();
#define TOCK()                                                            \
  {                                                                       \
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( \
        std::chrono::system_clock::now() - ____NOW____);                  \
    std::cout << elapsed.count() << "ms" << std::endl;                    \
  }

int main(int argc, char* argv[]) {
  CHECK(argc == 3, "Usage: ./a.out A_file.mtx B_file.mtx");

  std::string A_file(argv[1]), b_file(argv[2]);
  std::unique_ptr<CSSparseMatrix<double>> A =
      std::move(CSSparseMatrix<double>::ConstructFromFile(A_file));
  std::unique_ptr<DenseMatrix<double>> b = std::move(
      DenseMatrix<double>::ConstructFromFile(b_file, DenseStorageType::RMaj));

  {
    std::unique_ptr<DenseMatrix<double>> b_clone = std::move(b->clone());
    TICK();
    math::SimpleSparseTriangularSolve<double>(A.get(), b_clone.get());
    /*for (size_t i = 0; i < b_clone->rows(); ++i) {
      for (size_t j = 0; j < b_clone->cols(); ++j) {
        std::cout << *(b_clone->access(i, j)) << ' ';
      }
      std::cout << '\n';
    }*/
    TOCK();
  }

  {
    std::unique_ptr<DenseMatrix<double>> b_clone = std::move(b->clone());
    TICK();
    math::Opt1SparseTriangularSolve<double>(A.get(), b_clone.get());
    /*for (size_t i = 0; i < b_clone->rows(); ++i) {
      for (size_t j = 0; j < b_clone->cols(); ++j) {
        std::cout << *(b_clone->access(i, j)) << ' ';
      }
      std::cout << '\n';
    }*/
    TOCK();
  }
  return 0;
}
