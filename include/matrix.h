/*
 * Matrix data wrappers
 */
#include <cstring>
#include <map>
#include <memory>
#include <sstream>
#include <typeinfo>
#include <unordered_map>

#include "check_assert.h"
#include "mm_io.h"
#include "sp_elem_ptr.h"

enum SparseStorageType { CSC, CSR };
enum DenseStorageType { RMaj, CMaj };

template <class DType>
class Matrix {
 protected:
  std::unique_ptr<DType[]> data_;
  const size_t rows_, cols_;
  const size_t len_data_;

 public:
  Matrix(size_t rows, size_t cols, DType *data, size_t len_data)
      : rows_(rows), cols_(cols), data_(data), len_data_(len_data) {}

  size_t rows() const { return rows_; }
  size_t cols() const { return cols_; }
  DType *data() { return data_.get(); }
  size_t len_data() { return len_data; }

  virtual utils::SparseElementPointer<DType> access(size_t i, size_t j) = 0;
};

template <class DType>
class DenseMatrix : public Matrix<DType> {
 private:
  const DenseStorageType dense_type_;

  DenseMatrix(size_t rows, size_t cols, DType *data,
              DenseStorageType dense_type)
      : Matrix<DType>(rows, cols, data, rows * cols), dense_type_(dense_type) {}

 public:
  static std::unique_ptr<DenseMatrix> ConstructFromFile(
      std::string file_path, DenseStorageType dense_type) {
    if (typeid(double) != typeid(DType) && typeid(float) != typeid(DType)) {
      throw std::string("Only float and double elementary types supported");
    }

    const char *fname = file_path.c_str();
    int M, N;
    DType *val;
    bool r_maj = (dense_type == DenseStorageType::RMaj),
         is_double = (typeid(double) == typeid(DType));
    int ret_code = utils::mm_io::mm_read_array(fname, &M, &N, (char **)&val,
                                               is_double, r_maj);
    CHECK_EQ(ret_code, 0);

    size_t rows = M, cols = N;
    DType *data = new DType[rows * cols];
    for (size_t iter = 0; iter < (M * N); ++iter) {
      data[iter] = val[iter];
    }
    // memcpy(data, val, sizeof(DType) * rows * cols);
    DenseMatrix<DType> *ptr = new DenseMatrix(rows, cols, data, dense_type);
    std::unique_ptr<DenseMatrix<DType>> u_ptr(ptr);

    free(val);
    return std::move(u_ptr);
  }

  utils::SparseElementPointer<DType> access(size_t i, size_t j) {
    CHECK_LT(i, Matrix<DType>::rows_);
    CHECK_LT(j, Matrix<DType>::cols_);

    if (dense_type_ == DenseStorageType::RMaj) {
      return &(Matrix<DType>::data_[i * Matrix<DType>::cols_ + j]);
    } else {
      return &(Matrix<DType>::data_[i + j * Matrix<DType>::rows_]);
    }
  }

  inline DType *access_fast_RMaj(size_t i, size_t j) {
    return &(Matrix<DType>::data_[i * Matrix<DType>::cols_ + j]);
  }

  inline DType *access_fast_CMaj(size_t i, size_t j) {
    return &(Matrix<DType>::data_[i + j * Matrix<DType>::rows_]);
  }
};

template <class DType>
class SparseMatrix : public Matrix<DType> {
 protected:
  const SparseStorageType sparse_type_;
  size_t nnz_;

 public:
  SparseMatrix(size_t rows, size_t cols, DType *data, size_t len_data,
               SparseStorageType sparse_type, size_t nnz)
      : Matrix<DType>(rows, cols, data, len_data),
        sparse_type_(sparse_type),
        nnz_(nnz) {}
};

/*
 * Compressed Sparse storage specialization of SparseMatrix.
 *
 * NOTE: The sparsity of the matrix cannot increase in this implementation
 * NOTE: Storage currently limited to CSC format.
 */
template <class DType>
class CSSparseMatrix : public SparseMatrix<DType> {
 private:
  std::unique_ptr<size_t[]> I_data_, J_data_;

  CSSparseMatrix(size_t rows, size_t cols, DType *data,
                 SparseStorageType sparse_type, size_t nnz, size_t *I_data,
                 size_t *J_data)
      : SparseMatrix<DType>(rows, cols, data, nnz, sparse_type, nnz),
        I_data_(I_data),
        J_data_(J_data) {
    CHECK_EQ(sparse_type != SparseStorageType::CSC,
             "Currently handling only CSC format !");
  }

 public:
  static std::unique_ptr<CSSparseMatrix> ConstructFromFile(
      std::string file_path) {
    if (typeid(double) != typeid(DType) && typeid(float) != typeid(DType)) {
      throw std::string("Only float and double elementary types supported");
    }

    const char *fname = file_path.c_str();
    int M, N, nz;
    double *val;
    int *I, *J;

    int ret_code = utils::mm_io::mm_read_unsymmetric_sparse(fname, &M, &N, &nz,
                                                            &val, &I, &J);
    CHECK_EQ(ret_code, 0);

    std::unordered_map<size_t, std::unordered_map<size_t, double>>
        sparse_matrix_col_maj;
    for (int iter = 0; iter < nz; ++iter) {
      double value = val[iter];
      int i = I[iter];
      int j = J[iter];
      auto &col = sparse_matrix_col_maj[j];
      bool inserted = col.insert({i, value}).second;
      CHECK(inserted, "Multiple values for same location !");
    }

    free(val);
    free(I);
    free(J);

    std::map<size_t, std::map<size_t, double>> sparse_matrix_col_maj_ordered;
    for (auto &col_pair : sparse_matrix_col_maj) {
      size_t j = col_pair.first;
      auto &col = col_pair.second;
      auto &col_ordered = sparse_matrix_col_maj_ordered[j];
      for (auto &row_pair : col) {
        size_t i = row_pair.first;
        double value = row_pair.second;
        col_ordered[i] = value;
      }
    }

    size_t rows = M, cols = N, nnz = nz;
    DType *data = new DType[nnz];
    size_t *I_data = new size_t[nnz], *J_data = new size_t[cols + 1];

    for (size_t j = 0; j <= cols; ++j) {
      J_data[j] = 0;
    }
    size_t previous_col = 0;
    size_t data_iter = 0;
    for (auto &col_pair : sparse_matrix_col_maj_ordered) {
      size_t curr_col_num = col_pair.first;
      for (size_t j = previous_col + 2; j <= curr_col_num; ++j) {
        J_data[j] = J_data[j - 1];
      }
      auto &col = col_pair.second;
      for (auto &row_pair : col) {
        size_t curr_row_num = row_pair.first;
        double value = row_pair.second;
        I_data[data_iter] = curr_row_num;
        data[data_iter] = (DType)value;
        ++data_iter;
      }
      J_data[curr_col_num + 1] = J_data[curr_col_num] + col.size();
      previous_col = curr_col_num;
    }

    for (size_t j = previous_col + 2; j <= cols; ++j) {
      J_data[j] = J_data[j - 1];
    }

    auto *ptr = new CSSparseMatrix(rows, cols, data, SparseStorageType::CSC,
                                   nnz, I_data, J_data);
    std::unique_ptr<CSSparseMatrix> u_ptr(ptr);
    return std::move(u_ptr);
  }

  utils::SparseElementPointer<DType> access(size_t i, size_t j) {
    CHECK_LT(i, Matrix<DType>::rows_);
    CHECK_LT(j, Matrix<DType>::cols_);

    // NOTE: Assumes CSC storage format
    size_t col_idx = J_data_[j];
    size_t nxt_col_idx = J_data_[j + 1];

    // Currently using linear search since number of indices seem small
    for (size_t final_idx = col_idx; final_idx < nxt_col_idx; ++final_idx) {
      if (I_data_[final_idx] == i) {
        return &Matrix<DType>::data_[final_idx];
      }
      /*
       * if (I_data[final_idx] > i ) {
       *   return 0;
       * }
       */
    }

    // This is reached if,
    // 1. column j is empty, or
    // 2. ith row in jth column is empty
    return NULL;
  }

  size_t *I_data() const { return I_data_; }
  size_t *J_data() const { return I_data_; }

  CSSparseMatrix() = delete;
};
