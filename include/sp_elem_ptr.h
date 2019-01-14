/*
 * Wrapper to simplify access to sparse elements
 */
#ifndef _SP_ELEM_PTR_
#define _SP_ELEM_PTR_

namespace utils {

template <class DType>
class SparseElementPointer {
 private:
  DType* ptr_;

 public:
  SparseElementPointer(DType* ptr) : ptr_(ptr) {}
  DType operator*() const { return (ptr_ == NULL ? 0 : *ptr_); }
  DType* get() const { return ptr_; }
  bool is_null() { return ptr_ == nullptr; }

  SparseElementPointer() = delete;
};

}  // namespace utils

#endif
