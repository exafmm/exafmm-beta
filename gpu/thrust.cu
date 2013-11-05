#include <thrust/device_ptr.h>
#include <thrust/sort.h>

void sort(const int size, uint64_t * key, int * value) {
  thrust::device_ptr<uint64_t> keyBegin(key);
  thrust::device_ptr<uint64_t> keyEnd(key+size);
  thrust::device_ptr<int> valueBegin(value);
  thrust::sort_by_key(keyBegin, keyEnd, valueBegin);
}