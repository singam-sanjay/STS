#include <iostream>
#include <string>

#include "matrix.h"

int main(int argc, char* argv[]) {
  CHECK(argc == 3, "Need 2 arguments: ./a.out dense.mtx is_CMaj");

  std::string file_name(argv[1]);
  bool is_RMaj = (strcmp(argv[2], "0") == 0);
  bool is_CMaj = (strcmp(argv[2], "1") == 0);
  CHECK(is_RMaj || is_CMaj, "Invalid input for dense storage format");
  DenseStorageType storage_type =
      (is_RMaj ? DenseStorageType::RMaj : DenseStorageType::CMaj);

  std::unique_ptr<DenseMatrix<double>> mat = std::move(
      DenseMatrix<double>::ConstructFromFile(file_name, storage_type));

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
