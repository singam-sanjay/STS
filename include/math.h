#ifndef _MATH_
#define _MATH_

#include "check_assert.h"
#include "matrix.h"

namespace math {

/*
 * Unoptimized solver for triangular systems stored in a sparse matrix.
 *
 * A: lower triangular matrix stored in CSC format
 * b: column vector
 */
template <class DType>
void SimpleSparseTriangularSolve(CSSparseMatrix<DType>* A,
                                 DenseMatrix<DType>* b) {
  const size_t A_rows = A->rows(), A_cols = A->cols(), b_rows = b->rows(),
               b_cols = b->cols();
  CHECK_EQ(A_rows, b_rows);
  CHECK_EQ(b_cols, 1);

  DType* b_data = b->data();
  DType* A_data = A->data();
  size_t *I_data = A->I_data(), *J_data = A->J_data();
  for (size_t j = 0; j < A_cols; ++j) {
    b_data[j] /= A_data[J_data[j]];

    for (size_t i_idx = J_data[j] + 1; i_idx < J_data[j + 1]; ++i_idx) {
      REPORT(i_idx);
      b_data[I_data[i_idx]] -= A_data[i_idx] * b_data[j];
    }
  }
}

}  // namespace math

#endif
