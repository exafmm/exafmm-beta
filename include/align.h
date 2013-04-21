#ifndef align_h
#define align_h
#include <malloc.h>
#include <memory>

template <typename T, size_t NALIGN>
struct AlignedAllocator : public std::allocator<T> {
  using typename std::allocator<T>::size_type;
  using typename std::allocator<T>::pointer;

  template <typename U>
  struct rebind {
    typedef AlignedAllocator<U, NALIGN> other;
  };

  pointer allocate(size_type n) {
    void *ptr = NULL;
    int rc = posix_memalign(&ptr, NALIGN, n * sizeof(T));
    if (rc != 0) return NULL;
    if (ptr == NULL) throw std::bad_alloc();
    return reinterpret_cast<pointer>(ptr);
  }

  void deallocate(pointer p, size_type) {
    return free(p);
  }
};
#endif
