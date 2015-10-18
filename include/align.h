#ifndef align_h
#define align_h
#include <cstdlib>
#include <memory>

namespace exafmm {
  template <typename T, size_t nalign>
  struct AlignedAllocator : public std::allocator<T> {
    template <typename U>
    struct rebind {
      typedef AlignedAllocator<U, nalign> other;
    };

    T * allocate(size_t n) {
      void *ptr = NULL;
      int rc = posix_memalign(&ptr, nalign, n * sizeof(T));
      if (rc != 0) return NULL;
      if (ptr == NULL) throw std::bad_alloc();
      return reinterpret_cast<T*>(ptr);
    }

    void deallocate(T * p, size_t) {
      return free(p);
    }
  };
}
#endif
