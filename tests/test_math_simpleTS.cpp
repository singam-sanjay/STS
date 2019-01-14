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
#define TOCK_RETURN()                                      \
  ((std::chrono::duration_cast<std::chrono::milliseconds>( \
        std::chrono::system_clock::now() - ____NOW____))   \
       .count())

template <class DType>
void benchmark(char* argv[]) {
  std::string A_file(argv[2]), b_file(argv[3]);
  std::unique_ptr<CSSparseMatrix<DType>> A =
      std::move(CSSparseMatrix<DType>::ConstructFromFile(A_file));
  std::unique_ptr<DenseMatrix<DType>> b = std::move(
      DenseMatrix<DType>::ConstructFromFile(b_file, DenseStorageType::RMaj));

  std::unique_ptr<DenseMatrix<DType>> b_clone_1 = std::move(b->clone());
  {
    TICK();
    math::SimpleSparseTriangularSolve<DType>(A.get(), b_clone_1.get());
    TOCK();
    // b_clone_1->print(std::cout);
  }

  std::unique_ptr<DenseMatrix<DType>> b_clone_2 = std::move(b->clone());
  {
    TICK();
    math::Opt1SparseTriangularSolve<DType>(A.get(), b_clone_2.get());
    TOCK();
    // b_clone_2->print(std::cout);
  }

  // Verification
  std::cout << (*b_clone_1 == *b_clone_2 ? "Matrices equal !"
                                         : "Matrices not equal !")
            << std::endl;

  std::unique_ptr<DenseMatrix<DType>> b_clone_3 = std::move(b->clone());
  {
    TICK();
    math::Opt2SparseTriangularSolve<DType>(A.get(), b_clone_3.get());
    TOCK();
    // b_clone_2->print(std::cout);
  }

  // Verification
  std::cout << (*b_clone_1 == *b_clone_3 ? "Matrices equal !"
                                         : "Matrices not equal !")
            << std::endl;
  size_t time_for_5_iter = 0;
  for (size_t i = 0; i < 15; ++i) {
    std::unique_ptr<DenseMatrix<DType>> b_local_clone = std::move(b->clone());
    TICK();
    math::SimpleSparseTriangularSolve<DType>(A.get(), b_local_clone.get());
    if (i > 9) {
      time_for_5_iter += TOCK_RETURN();
    }
  }
  std::cout << "Average time (naive) : " << time_for_5_iter / 5.0 << '\n';

  time_for_5_iter = 0;
  for (size_t i = 0; i < 15; ++i) {
    std::unique_ptr<DenseMatrix<DType>> b_local_clone = std::move(b->clone());
    TICK();
    math::Opt1SparseTriangularSolve<DType>(A.get(), b_local_clone.get());
    if (i > 9) {
      time_for_5_iter += TOCK_RETURN();
    }
  }
  std::cout << "Average time (Opt1) : " << time_for_5_iter / 5.0 << '\n';

  time_for_5_iter = 0;
  for (size_t i = 0; i < 15; ++i) {
    std::unique_ptr<DenseMatrix<DType>> b_local_clone = std::move(b->clone());
    TICK();
    math::Opt2SparseTriangularSolve<DType>(A.get(), b_local_clone.get());
    if (i > 9) {
      time_for_5_iter += TOCK_RETURN();
    }
  }
  std::cout << "Average time (Opt2) : " << time_for_5_iter / 5.0 << '\n';
}

int main(int argc, char* argv[]) {
  CHECK(argc == 4, "Usage: ./a.out is_double A_file.mtx B_file.mtx");

  if (strcmp(argv[1], "1") == 0) {
    benchmark<double>(argv);
  } else if (strcmp(argv[1], "0") == 0) {
    benchmark<float>(argv);
  } else {
    CHECK(false, "Invalid option for is_double.");
  }
  return 0;
}
