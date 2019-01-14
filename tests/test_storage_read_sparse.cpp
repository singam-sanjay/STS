#include <iostream>
#include <string>

#include "matrix.h"

int main(int argc, char* argv[]) {
  CHECK(argc != 2, "Need 1 argument: ./a.out sparse.mtx");

  std::string file_name(argv[1]);
  std::unique_ptr<CSSparseMatrix<double>> mat =
      std::move(CSSparseMatrix<double>::ConstructFromFile(file_name));

  for (size_t i = 0; i < mat->rows(); ++i) {
    for (size_t j = 0; j < mat->cols(); ++j) {
      auto ptr = mat->access(i, j);
      std::cout << (ptr.is_null() ? std::string("NULL") : std::to_string(*ptr))
                << ' ';
    }
    std::cout << std::endl;
  }
  return 0;
}
