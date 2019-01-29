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
      b_data[I_data[i_idx]] -= A_data[i_idx] * b_data[j];
    }
  }
}

/*
 * Level 1 optimized solver for triangular systems stored in a sparse matrix.
 *
 * A: lower triangular matrix stored in CSC format
 * b: column vector
 */
template <class DType>
void Opt1SparseTriangularSolve(CSSparseMatrix<DType>* A,
                               DenseMatrix<DType>* b) {
  const size_t A_rows = A->rows(), A_cols = A->cols(), b_rows = b->rows(),
               b_cols = b->cols();
  CHECK_EQ(A_rows, b_rows);
  CHECK_EQ(b_cols, 1);

  DType* b_data = b->data();
  DType* A_data = A->data();
  size_t *I_data = A->I_data(), *J_data = A->J_data();

#define VECTOR_SIZE 8
  DType vector_b[VECTOR_SIZE], vector_A[VECTOR_SIZE];
  int indices[VECTOR_SIZE];
  for (size_t j = 0; j < A_cols; ++j) {
    DType b_solution = b_data[j];
    const size_t curr_col_idx = J_data[j], next_col_idx = J_data[j + 1];
    b_solution /= A_data[curr_col_idx];
    b_data[j] = b_solution;

    /*
    // Not vectorized: according to -fopt-info-vec-all
    #pragma unroll
    for (size_t i_idx = curr_col_idx + 1; i_idx < next_col_idx; ++i_idx) {
      b_data[I_data[i_idx]] -= A_data[i_idx] * b_solution;
    }
    */

    const int upper_bound = ((int)next_col_idx - 1);
    const int lower_bound = ((int)curr_col_idx + 1);
    const int num_iter = upper_bound - lower_bound + 1;
    int i;
    for (i = 0; (i + 7) < num_iter; i += VECTOR_SIZE) {
      for (int ii = 0; ii < VECTOR_SIZE; ++ii) {
        indices[ii] = I_data[lower_bound + i + ii];
        vector_b[ii] = b_data[indices[ii]];
        vector_A[ii] = A_data[lower_bound + i + ii];
      }

      // NOTE: GCC doesn't vectorize this !
      for (int ii = 0; ii < VECTOR_SIZE; ++ii) {
        vector_b[ii] -= (vector_A[ii] * b_solution);
      }

      for (int ii = 0; ii < VECTOR_SIZE; ++ii) {
        b_data[indices[ii]] = vector_b[ii];
      }
    }

    for (/*i = i*/; i < num_iter; ++i) {
      int index = I_data[lower_bound + i];
      DType b_elem = b_data[index];
      b_elem -= A_data[lower_bound + i] * b_solution;
      b_data[index] = b_elem;
    }
  }
}

/*
 * Level 2 optimized solver for triangular systems stored in a sparse matrix.
 *
 * A: lower triangular matrix stored in CSC format
 * b: column vector
 */
template <class DType>
void Opt2SparseTriangularSolve(CSSparseMatrix<DType>* A,
                               DenseMatrix<DType>* b) {
  const size_t A_rows = A->rows(), A_cols = A->cols(), b_rows = b->rows(),
               b_cols = b->cols();
  CHECK_EQ(A_rows, b_rows);
  CHECK_EQ(b_cols, 1);

  DType* b_data = b->data();
  DType* A_data = A->data();
  size_t *I_data = A->I_data(), *J_data = A->J_data();

#pragma omp parallel num_threads(2)
  {
#define VECTOR_SIZE 8
    // DType vector_b[VECTOR_SIZE], vector_A[VECTOR_SIZE];
    // int indices[VECTOR_SIZE];
    for (size_t j = 0; j < A_cols; ++j) {
      double b_solution = b_data[j];
      const size_t curr_col_idx = J_data[j], next_col_idx = J_data[j + 1];
      b_solution /= A_data[curr_col_idx];
#pragma omp barrier
#pragma omp single
      { b_data[j] = b_solution; }

      /*
      // Not vectorized: according to -fopt-info-vec-all
      #pragma unroll
      for (size_t i_idx = curr_col_idx + 1; i_idx < next_col_idx; ++i_idx) {
        b_data[I_data[i_idx]] -= A_data[i_idx] * b_solution;
      }
      */

      const int upper_bound = ((int)next_col_idx - 1);
      const int lower_bound = ((int)curr_col_idx + 1);
      const int num_iter = upper_bound - lower_bound + 1;
      /*
      int i;
      #pragma omp parallel for nowait
      for (i = 0 ; (i+7) < num_iter ; i += VECTOR_SIZE) {
        for( int ii = 0 ; ii < VECTOR_SIZE ; ++ii ) {
          indices[ii] = I_data[lower_bound+i+ii];
          vector_b[ii] = b_data[indices[ii]];
          vector_A[ii] = A_data[lower_bound+i+ii];
        }

        for( int ii = 0 ; ii < VECTOR_SIZE ; ++ii ) {
          vector_b[ii] -= (vector_A[ii] * b_solution);
        }

        for( int ii = 0 ; ii < VECTOR_SIZE ; ++ii ) {
          b_data[indices[ii]] = vector_b[ii];
        }
      }
      */

#pragma omp for
      for (int i = 0; i < num_iter; ++i) {
        int index = I_data[lower_bound + i];
        DType b_elem = b_data[index];
        b_elem -= A_data[lower_bound + i] * b_solution;
        b_data[index] = b_elem;
      }
    }
  }
}

}  // namespace math

#endif
