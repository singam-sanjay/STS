/*
 * Matrix data wrappers
 */
#include<memory>
#include<sstream>

enum SparseStorageType { CSC, CSR };

template <class DType>
class Matrix {
 protected:
  std::unique_ptr<DType[]> data_;
  const size_t rows_, cols_;
 
 public:
  Matrix(size_t rows, size_t cols, DType *data)
      : rows_(rows), cols_(cols), data_(data) {}
  
  size_t rows() const { return rows_; }
  size_t cols() const { return cols_; }
  DType* data() { return data_.get(); }

  virtual DType& access( size_t i, size_t j);
};

template <class DType>
class SparseMatrix : public Matrix<DType> {
 protected:
  const SparseStorageType sparse_type_;
  size_t nnz_;

 public:
  SparseMatrix(size_t rows, size_t cols, DType *data,
               SparseStorageType sparse_type, size_t nnz)
      : Matrix<DType>(rows, cols, data), sparse_type_(sparse_type), nnz_(nnz) {}
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
 
 public:
  CSSparseMatrix(size_t rows, size_t cols, DType *data,
                 SparseStorageType sparse_type, size_t nnz)
      : SparseMatrix<DType>(rows, cols, data, sparse_type, nnz) {
    if (sparse_type != SparseStorageType::CSC) {
      throw std::string("Currently handling only CSC format !");
    }
  }

  static std::unique_ptr<DType> ConstructFromFile(std::string file_path);

  DType& access(size_t i, size_t j) {
    // NOTE: Assumes CSC storage format
    if( i > Matrix<DType>::rows_ ) {
      std::stringstream ss;
      ss << "Row #" << i << "requested in matrix with " << Matrix<DType>::rows_ << " rows.";
      throw ss.str();
    }
    if( j > Matrix<DType>::cols_ ) {
      std::stringstream ss;
      ss << "Col #" << j << "requested in matrix with " << Matrix<DType>::cols_ << " cols.";
      throw ss.str();
    }

    size_t col_idx = J_data_[j];
    size_t nxt_col_idx = J_data_[j+1];

    // Currently using linear search since number of indices seem small
    for( size_t final_idx = col_idx ; final_idx < nxt_col_idx ; ++final_idx ) {
      if( I_data_[final_idx] == i ) {
        return Matrix<DType>::data_[final_idx];
      }
    }

    return 0;
  }
};
